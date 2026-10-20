%% Test open-loop for the quadcopter in a Spiral Trajectory — B01 model.
%
%  Spiral (elliptical helix, climbing) trajectory tracking using
%  differential flatness.
%
%  Frames:
%       Inertial : NED
%       Body     : FRD
%       Velocity state : NED (v_n, v_e, v_d).

clear; clc; close all;

%% PATHS
addpath(fullfile('models'));
addpath(fullfile('flatness'));
addpath(fullfile('configs'));
addpath(fullfile('trajectories'));
addpath(fullfile('utils'));

%% MODEL SELECTION
mp = ModelParameters;
tp = TrajectoryParameters;

quad = B01_ANTX_quadcopter(mp, false);   % delay_motors 
flatness_map = B01_FlatnessMap;
fprintf('Selected model: B01 (16-State, rotational damping + motor delay)\n');

fprintf('OpenLoopSpiral \n');
fprintf('Drone: m=%.3f kg, T_hover_total=%.3f N, T_hover_per_motor=%.4f N\n\n', ...
    mp.m, mp.m*mp.g, mp.m*mp.g/4);

%% TRAJECTORY CONFIGURATION
sp = tp.Spiral;

traj = struct( ...
    'name', 'Spiral', ...
    'obj',  ShapeSpiral(sp), ...
    'col',  [0.9 0.55 0.1]);

T_final = 1 * tp.T_duration_spiral;
fprintf('Simulation time: %.3f s\n\n', T_final);

%% ODE OPTIONS
opts_ode = odeset('RelTol', 1e-7, 'AbsTol', 1e-8);

%% INITIAL CONDITION FROM FLATNESS MAP
n_states = 12 + 4*quad.delay_motors;

[s0, ds0, dds0, ddds0, dddds0] = traj.obj.get_flat_outputs(0);
[x0_full, u0] = flatness_map.map(mp, s0, ds0, dds0, ddds0, dddds0);
x0 = x0_full(1:n_states);

init_phi = x0(7); init_theta = x0(8); init_psi = x0(9);
init_p   = x0(10); init_q = x0(11); init_r = x0(12);

fprintf('Initial condition:\n');
fprintf('  x_I   = %.2f m\n', x0(1));
fprintf('  y_I   = %.2f m\n', x0(2));
fprintf('  z_I   = %.2f m\n', x0(3));
fprintf('  v_n   = %.2f m/s\n', x0(4));
fprintf('  v_e   = %.2f m/s\n', x0(5));
fprintf('  v_d   = %.2f m/s\n', x0(6));
fprintf('  phi   = %.2f deg\n', init_phi*180/pi);
fprintf('  theta = %.2f deg\n', init_theta*180/pi);
fprintf('  psi   = %.2f deg\n', init_psi*180/pi);
fprintf('  p     = %.2f rad/s\n', init_p);
fprintf('  q     = %.2f rad/s\n', init_q);
fprintf('  r     = %.2f rad/s\n\n', init_r);

fprintf('Initial motor thrusts (x0):\n');
if quad.delay_motors
    fprintf('  T1 = %.4f N   T2 = %.4f N   T3 = %.4f N   T4 = %.4f N\n\n', ...
        x0(13), x0(14), x0(15), x0(16));
else
    fprintf('  (no motor-delay states in this mode)\n\n');
end

fprintf('Initial inputs (desired motor thrusts):\n');
fprintf('  T1_d = %.4f N   T2_d = %.4f N   T3_d = %.4f N   T4_d = %.4f N\n\n', ...
    u0(1), u0(2), u0(3), u0(4));

%% ODE INTEGRATION
[t_sol, x_sol] = ode45( ...
    @(t, x) ode_dynamics(t, x, quad, traj.obj, mp, flatness_map), ...
    [0 T_final], x0, opts_ode);

fprintf('Final integration time: %.4f s\n\n', t_sol(end));

%% RECONSTRUCT REFERENCES
N = length(t_sol);
x_ref = zeros(N, n_states);
u_ref = zeros(N, 4);

for k = 1:N
    [s, ds, dds, ddds, dddds] = traj.obj.get_flat_outputs(t_sol(k));
    [xr_full, ur] = flatness_map.map(mp, s, ds, dds, ddds, dddds);
    x_ref(k, :) = xr_full(1:n_states)';
    u_ref(k, :) = ur';
end

% Unwrap yaw angle to handle rotation discontinuities
x_ref(:, 9) = unwrap(x_ref(:, 9));
x_sol(:, 9) = unwrap(x_sol(:, 9));

%% POSITION ERROR
pos_err = vecnorm(x_sol(:, 1:3) - x_ref(:, 1:3), 2, 2);
fprintf('Max position error: %.3e m\n', max(pos_err));

% Motor thrust tracking error
if quad.delay_motors
    thrust_err = vecnorm(x_sol(:, 13:16) - x_ref(:, 13:16), 2, 2);
    fprintf('Max motor-thrust-state error: %.3e N\n', max(thrust_err));
end

%% FIGURE 1 — RIGID-BODY STATE TRACKING (states 1-12)
state_names = {
    'r_n [m]','r_e [m]','r_d [m]',...
    'v_n [m/s]','v_e [m/s]','v_d [m/s]',...
    '\phi [deg]','\theta [deg]','\psi [deg]',...
    'p [rad/s]','q [rad/s]','r [rad/s]'};

angle_scale = [1 1 1 1 1 1 180/pi 180/pi 180/pi 1 1 1];

figure('Name', 'Spiral State Tracking', 'Color', 'w', 'Position', [30 50 1600 850]);
tiledlayout(4, 3, 'Padding', 'compact', 'TileSpacing', 'compact');

for j = 1:12
    nexttile;
    plot(t_sol, x_ref(:, j)*angle_scale(j), 'r--', 'LineWidth', 1.8); hold on;
    plot(t_sol, x_sol(:, j)*angle_scale(j), 'b', 'LineWidth', 1.2);
    ylabel(state_names{j}); grid on;

    if j <= 3, title(state_names{j}); end
    if j > 9,  xlabel('t [s]'); end
    if j == 1, legend('ref', 'sim', 'Location', 'best'); end
end
sgtitle(sprintf('Spiral tracking (rigid body)  |  max pos error = %.2e m', max(pos_err)));

%% FIGURE 1b — MOTOR THRUST STATES (13-16)
T_hover_pm = mp.m * mp.g / 4;
if quad.delay_motors
    figure('Name', 'Spiral Motor Thrusts', 'Color', 'w', 'Position', [30 950 1200 350]);
    motor_labels = {'T_1 [N]', 'T_2 [N]', 'T_3 [N]', 'T_4 [N]'};
    for j = 1:4
        subplot(1, 4, j);
        plot(t_sol, x_ref(:, 12+j), 'r--', 'LineWidth', 1.8); hold on;
        plot(t_sol, x_sol(:, 12+j), 'b', 'LineWidth', 1.2);
        yline(T_hover_pm, 'k:');
        ylabel(motor_labels{j}); xlabel('t [s]'); grid on;
        if j == 1, legend('ref (no delay)', 'sim (with delay)', 'Location', 'best'); end
    end
    sgtitle(sprintf('Motor thrust tracking  |  max err = %.2e N', max(thrust_err)));
end

%% FIGURE 2 — INPUTS (desired motor thrusts)
figure('Name', 'Spiral Inputs', 'Color', 'w', 'Position', [30 50 900 500]);
input_labels = {'T_1^d [N]', 'T_2^d [N]', 'T_3^d [N]', 'T_4^d [N]'};

for j = 1:4
    subplot(4, 1, j);
    plot(t_sol, u_ref(:, j), 'Color', traj.col, 'LineWidth', 1.8); grid on;
    yline(T_hover_pm, 'k--', 'T_{hover}/4');
    ylabel(input_labels{j});
end
xlabel('t [s]'); sgtitle('Reference inputs (desired motor thrusts)');

%% FIGURE 3 — POSITION ERROR
figure('Name', 'Position Error', 'Color', 'w');
plot(t_sol, pos_err, 'LineWidth', 2); grid on;
xlabel('t [s]'); ylabel('||e_p|| [m]');
title('Position tracking error');

%% FIGURE 4 — 3D DRONE ANIMATION
fig3 = figure('Name', '3D Spiral Animation', 'Color', 'w', 'Position', [900 100 800 750]);
ax3 = axes(fig3);
PlotDrone2.setup_axes(ax3, '3D Drone — Spiral Trajectory (B01)');
hold(ax3, 'on');

sz = 2.5;
PlotDrone2.draw_ground(ax3, [-sz sz], [-sz sz], 0);

plot3(ax3, x_ref(:, 1), x_ref(:, 2), x_ref(:, 3), 'k--', 'LineWidth', 1.5);
h_traj = plot3(ax3, NaN, NaN, NaN, 'Color', traj.col, 'LineWidth', 1.5);
h_drone = PlotDrone2.init(ax3, mp);

pos0 = x_sol(1, 1:3)';
R0 = Euler2R_ZYX(x_sol(1, 7), x_sol(1, 8), x_sol(1, 9));
PlotDrone2.update(h_drone, pos0, R0);

x_all = [x_ref(:, 1); x_sol(:, 1)];
y_all = [x_ref(:, 2); x_sol(:, 2)];
z_all = [x_ref(:, 3); x_sol(:, 3)];
margin = 0.5;

xlim(ax3, [min(x_all)-margin max(x_all)+margin]);
ylim(ax3, [min(y_all)-margin max(y_all)+margin]);
zlim(ax3, [min(z_all)-margin max(z_all)+margin]);

annotation(fig3, 'textbox', [0.01 0.01 0.28 0.14],...
    'String', { ...
    'FRD body frame:', ...
    ' x_B (red)   = Forward', ...
    ' y_B (green) = Right', ...
    ' z_B (blue)  = Down', ...
    ' \phi > 0 : right wing DOWN', ...
    ' \theta > 0 : nose UP', ...
    ' \psi > 0 : yaw right'}, ...
    'FontSize', 7, 'BackgroundColor', 'w', 'EdgeColor', [0.5 0.5 0.5]);

skip = max(1, floor(length(t_sol)/140));

for k = 1:skip:length(t_sol)
    pos_k = x_sol(k, 1:3)';
    phi_k = x_sol(k, 7);
    the_k = x_sol(k, 8);
    psi_k = x_sol(k, 9);
    R_k   = Euler2R_ZYX(phi_k, the_k, psi_k);

    PlotDrone2.update(h_drone, pos_k, R_k);
    set(h_traj, 'XData', x_sol(1:k, 1), 'YData', x_sol(1:k, 2), 'ZData', x_sol(1:k, 3));

    title(ax3, sprintf([ ...
        'Spiral  |  t = %.2f s  |  ', ...
        '\\phi = %.1f°  ', ...
        '\\theta = %.1f°  ', ...
        '\\psi = %.1f°'], ...
        t_sol(k), phi_k*180/pi, the_k*180/pi, psi_k*180/pi));
    drawnow;
end

fprintf('\n End simulation \n');

%% LOCAL FUNCTIONS
function dxdt = ode_dynamics(t, x, quad, traj_obj, mp, flatness_map)
    % Reconstruct trajectory state and derivatives
    [s, ds, dds, ddds, dddds] = traj_obj.get_flat_outputs(t);

    % Get reference inputs (desired motor thrusts) from Flatness Map
    [~, uref] = flatness_map.map(mp, s, ds, dds, ddds, dddds);

    % Evaluate dynamics (16-state, motor delay included in quad.dynamics)
    dxdt = quad.dynamics(t, x, uref);
end