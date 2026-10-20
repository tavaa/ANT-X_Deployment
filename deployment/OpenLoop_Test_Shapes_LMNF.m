%% OpenLoop_Test_Shapes_LMNF - test open loop (una shape alla volta) sul drone ANT-X  DroneCmd
%
%  Shape: Circle (Deployment) | Roll | Yaw | Thrust (una per run, vedi SHAPE)
%  Pipeline: Shape -> FlatnessMap -> u_ref(t_k) =  [T L M N]
% -> normalizzazione LMNF -> sendLMNFsetpointOL
%
%  Frames:
%       Inertial : NED
%       Body     : FRD
%       Velocity state : NED (v_x, v_y, v_z).
%
%  Modello A02.
%  Plant: drone reale (Deployment).

clear;
clc;
close all;

%% PATHS
addpath(fullfile('models'));
addpath(fullfile('flatness'));
addpath(fullfile('configs'));
addpath(fullfile('trajectories'));
addpath(fullfile('utils'));
addpath(fullfile('deployment'));

%% PARAMETRI
mp = ModelParameters;
tp = TrajectoryParameters;

ns       = '/antx_dronelab';     % namespace ROS del drone (da adattare)
Ts       = 0.1;                  % [s] passo di campionamento = durata di ogni comando
max_tilt = deg2rad(20);          % [rad] soglia di sicurezza in open loop
p0       = [0; 0; -1];           % [m] setpoint di partenza (NED), uguale per tutte le shape
psi0     = 0;                    % [rad] heading iniziale

%% Traiettoria da testare
SHAPE = 'Circle';                % 'Circle' | 'Roll' | 'Yaw' | 'Thrust'
 
switch SHAPE
    case 'Circle',     shape = ShapeCircleDeployment(tp.CircleDeployment);  T_duration = shape.duration();
    case 'Roll',       shape = ShapeRollTest(tp.Roll);                      T_duration = tp.T_duration_roll;
    case 'Yaw',        shape = ShapeYawTest(tp.Yaw);                        T_duration = tp.T_duration_yaw;
    case 'Thrust',     shape = ShapeThrustTest(tp.Thrust);                  T_duration = tp.T_duration_thrust;
    otherwise, error('OpenLoop_Test_Shapes_LMNF:Shape', 'Unknown shape "%s".', SHAPE);
end

t_k        = 0:Ts:T_duration;
N          = numel(t_k) - 1;     % numero di comandi (ognuno tenuto per Ts)

%% REFERENCE: campionamento della mappa di flatness u_ref = u(t_k)
% Il riferimento e' algebrico.
% x_ref(:,1) = hover nel setpoint.

x_ref = zeros(12, N+1);
u_ref = zeros(4,  N+1);
for k = 1:N+1
    [s, ds, dds, ddds, dddds] = shape.get_flat_outputs(t_k(k));
    [x_ref(:,k), u_ref(:,k)]  = A02_FlatnessMap.map(mp, s, ds, dds, ddds, dddds);
end
x_ref(9,:) = unwrap(x_ref(9,:));  % la mappa da' psi in [-pi, pi]: unwrap lo rende continuo (0 -> 2*pi per giro)

%% NORMALIZZAZIONE LMNF
d    = mp.b/sqrt(2);              % braccio proiettato
Tmax = mp.K_T*mp.Omega_max^2;     % spinta massima del singolo rotore

% scale = massimo di ogni canale con i rotori in [0, Tmax].
LMNF = [ u_ref(2,:) / (2*d*Tmax);                      % L in [-1, 1]
         u_ref(3,:) / (2*d*Tmax);                      % M in [-1, 1]
         u_ref(4,:) / (2*mp.Beta*Tmax);                % N in [-1, 1]
         zeros(1, N+1);                            % Fx (non usato)
         zeros(1, N+1);                            % Fy (non usato)
         u_ref(1,:) / (4*Tmax) ];                      % Fz in [0, 1]

Fz_atteso = mp.m*mp.g / (4*Tmax);
fprintf('Durata %.1f s, %d comandi | Fz in [%.4f %.4f] | |L|,|M| max %.1e | |N| max %.1e\n', ...
    T_duration, N, min(LMNF(6,:)), max(LMNF(6,:)), max(max(abs(LMNF(1:2,:)))), max(abs(LMNF(3,:))));

%% SETUP DRONE: connessione, arm, decollo, posizionamento sul setpoint
% takeoff e sendPositionSetpoint ritornano subito: le pause servono ad aspettare
% che il drone sia davvero fermo in hover nel setpoint (= x_ref(:,1)).
Drone.startup();
drone = Drone(ns);
drone.connect();
drone.arm();
pause(2);

drone.takeoff(p0(3));
pause(5);
drone.sendPositionSetpoint(p0, psi0);
pause(5);                                          % attesa drone fermo in hover

%% MISURA FORZA VERTICALE IN HOVER (verifica delle scale)
[~, ~, Fz_hover] = drone.getAppliedForce();
fprintf('Fz in hover: misurato %.3f, atteso %.3f\n', Fz_hover, Fz_atteso);


%% MAIN LOOP: ogni Ts invia il comando k e registra lo stato
X_dr = nan(12, N);                 % stato misurato
t_dr = nan(1,  N);                 % istante della misura

t0 = tic;
for k = 1:N

    % comando k = u_ref(t_k); resta attivo fino al comando successivo (ZOH)
    drone.sendLMNFsetpointOL(LMNF(:,k)', true, true, true, false, false, true);

    % stato misurato subito dopo l'invio: corrisponde a x(t_k), confrontabile con x_ref(:,k)
    [pn, pe, pd] = drone.getPosition();
    [vn, ve, vd] = drone.getVelocity();
    [ph, th, ps] = drone.getEuler();
    [p,  q,  r ] = drone.getAngularVelocity();
    
    % costruisci vettore stato
    X_dr(:,k) = [pn; pe; pd; vn; ve; vd; ph; th; ps; p; q; r];
    t_dr(k)   = toc(t0);

    % sicurezza: in open loop l'assetto diverge, si interrompe prima del ribaltamento
    if acos(cos(ph)*cos(th)) > max_tilt
        fprintf('STOP: tilt oltre soglia a t = %.2f s\n', t_dr(k));
        break;
    end

    while toc(t0) < k*Ts, end                      % passo Ts in tempo reale
end
k_end = k;

%% FINE: stop iniezione, hover, atterraggio
drone.sendLMNFsetpointOL(zeros(1,6), false, false, false, false, false, false);  % torna il controllore di bordo

[pn, pe, pd] = drone.getPosition();
drone.sendPositionSetpoint([pn; pe; pd], psi0); % hover dove si trova
pause(3);
drone.land();
pause(5);                                          % disarm solo a terra
drone.disarm();
drone.disconnect();
Drone.shutdown();

%% RISULTATI: riferimento vs drone
X_dr      = X_dr(:, 1:k_end);
t_dr      = t_dr(1:k_end);
X_dr(9,:) = unwrap(X_dr(9,:));
e_p       = vecnorm(X_dr(1:3,:) - x_ref(1:3,1:k_end));

figure('Name', ['OpenLoop ' SHAPE ': riferimento vs drone'], 'Color', 'w');

subplot(2,2,[1 3]);
plot(x_ref(2,:), x_ref(1,:), 'k--', X_dr(2,:), X_dr(1,:), 'r', 'LineWidth', 1.3); hold on;
plot(p0(2), p0(1), 'ko', 'MarkerFaceColor', 'k');
axis equal; grid on; xlabel('E [m]'); ylabel('N [m]');
legend('riferimento', 'drone', 'setpoint'); title('Traiettoria XY');

subplot(2,2,2);
plot(t_dr, e_p, 'r', 'LineWidth', 1.3); grid on;
xlabel('t [s]'); ylabel('||e_p|| [m]'); title('Errore di posizione');

subplot(2,2,4);
plot(t_k, rad2deg(x_ref(7:8,:)), 'k--', t_dr, rad2deg(X_dr(7:8,:)), 'LineWidth', 1.3); grid on;
xlabel('t [s]'); ylabel('[deg]'); legend('\phi_{ref}', '\theta_{ref}', '\phi', '\theta');
title('Assetto');

sgtitle(sprintf(['%s | Open loop %.1f / %.1f s (%d / %d comandi) | Fz hover %.3f (atteso %.3f) | ' ...
                 'e_p max %.2f m'], SHAPE, t_dr(end), T_duration, k_end, N, Fz_hover, Fz_atteso, max(e_p)));

save(['OpenLoop_' SHAPE '_log.mat'], 't_k', 'x_ref', 'u_ref', 'LMNF', 't_dr', 'X_dr', 'p0', 'psi0', 'Fz_hover');