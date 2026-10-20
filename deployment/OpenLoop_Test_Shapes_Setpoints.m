%% OpenLoop_Shape - test open loop del cerchio sul drone ANT-X  DroneCmd
%
%  Pipeline: ShapeCircleDeployment -> FlatnessMap -> u_ref(t_k) =  [T L M N]
%  -> ode45 con A02 -> sendPositionSetpoint
%
%  Frames:
%       Inertial : NED
%       Body     : FRD
%       Velocity state : NED (v_x, v_y, v_z).
%
%  Modello A02.
%  Plant: drone reale (Deployment), guidato con setpoint di posizione:
%  ode45 + modello A02 propagano lo stato con u_ref(t_k) (ZOH, Ts) e la
%  posizione calcolata e' il setpoint del drone.

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
model    = A02_SimplifiedModel_Deployment(mp);   % modello A02 per ode45
ode_opt  = odeset('RelTol', 1e-10, 'AbsTol', 1e-12);
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

%% SETUP DRONE: connessione, arm, decollo, posizionamento sul setpoint
% takeoff e sendPositionSetpoint ritornano subito: le pause servono ad aspettare
% che il drone sia davvero fermo in hover nel setpoint (= x_ref(:,1)).
Drone.startup();
drone = Drone(ns);
drone.connect();
drone.arm();
pause(2);

drone.takeoff(cp.p0(3));
pause(5);
drone.sendPositionSetpoint(cp.p0, cp.psi0);
pause(5);                                          % attesa drone fermo in hover


%% MAIN LOOP: ogni Ts invia il comando k e registra lo stato
X_dr = nan(12, N);                 % stato misurato
t_dr = nan(1,  N);                 % istante della misura

x_sim = x_ref(1:12,1);             % stato del modello A02 (hover nel setpoint)

t0 = tic;
for k = 1:N

    % setpoint k = posizione del modello A02 a t_k
    drone.sendPositionSetpoint(x_sim(1:3), x_sim(9));

    % stato misurato subito dopo l'invio: corrisponde a x(t_k), confrontabile con x_ref(:,k)
    [pn, pe, pd] = drone.getPosition();
    [vn, ve, vd] = drone.getVelocity();
    [ph, th, ps] = drone.getEuler();
    [p,  q,  r ] = drone.getAngularVelocity();
    
    % costruisci vettore stato
    X_dr(:,k) = [pn; pe; pd; vn; ve; vd; ph; th; ps; p; q; r];
    t_dr(k)   = toc(t0);

    % prossimo stato: ode45 con A02, u_ref(t_k) tenuto costante per Ts (ZOH)
    [~, x_ode] = ode45(@(t, x) model.dynamics(t, x, u_ref(:,k)), [t_k(k), t_k(k)+Ts], x_sim, ode_opt);
    x_sim = x_ode(end,:)';

    % sicurezza: in open loop l'assetto diverge, si interrompe prima del ribaltamento
    if acos(cos(ph)*cos(th)) > max_tilt
        fprintf('STOP: tilt oltre soglia a t = %.2f s\n', t_dr(k));
        break;
    end

    while toc(t0) < k*Ts, end                      % passo Ts in tempo reale
end
k_end = k;

%% FINE: stop iniezione, hover, atterraggio
[pn, pe, pd] = drone.getPosition();
drone.sendPositionSetpoint([pn; pe; pd], cp.psi0); % hover dove si trova
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

figure('Name', 'OpenLoop Circle: riferimento vs drone', 'Color', 'w');

subplot(2,2,[1 3]);
plot(x_ref(2,:), x_ref(1,:), 'k--', X_dr(2,:), X_dr(1,:), 'r', 'LineWidth', 1.3); hold on;
plot(cp.p0(2), cp.p0(1), 'ko', 'MarkerFaceColor', 'k');
axis equal; grid on; xlabel('E [m]'); ylabel('N [m]');
legend('riferimento', 'drone', 'setpoint'); title('Traiettoria XY');

subplot(2,2,2);
plot(t_dr, e_p, 'r', 'LineWidth', 1.3); grid on;
xlabel('t [s]'); ylabel('||e_p|| [m]'); title('Errore di posizione');

subplot(2,2,4);
plot(t_k, rad2deg(x_ref(7:8,:)), 'k--', t_dr, rad2deg(X_dr(7:8,:)), 'LineWidth', 1.3); grid on;
xlabel('t [s]'); ylabel('[deg]'); legend('\phi_{ref}', '\theta_{ref}', '\phi', '\theta');
title('Assetto');

sgtitle(sprintf('Open loop %.1f / %.1f s (%d / %d comandi) | e_p max %.2f m', ...
    t_dr(end), T_duration, k_end, N, max(e_p)));

save('OpenLoop_Shape_log.mat', 't_k', 'x_ref', 'u_ref', 't_dr', 'X_dr', 'cp');