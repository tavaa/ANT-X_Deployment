%% Test open-loop for the quadcopter in a Lemniscate Trajectory.
%
%  Lemniscate (figure-8) trajectory tracking using differential flatness.
%
%  Expected:
%       - variable tilt (roll/pitch change continuously to track the 8)
%       - constant altitude
%       - yaw rotating to stay tangent to the trajectory
%
%  Frames:
%       Inertial : NED
%       Body     : FRD

clear; clc; close all;

%% PATHS
addpath(fullfile('models'));
addpath(fullfile('flatness'));
addpath(fullfile('configs'));
addpath(fullfile('trajectories'));
addpath(fullfile('utils'));

%% SETTINGS
use_quaternions = false; % Set to true for A01 (quaternion), false for A00 (Euler Z-Y-X)

%% DYNAMIC MODEL SELECTION
mp   = ModelParameters;
tp   = TrajectoryParameters;

if use_quaternions
    quad = A01_SimplifiedModel_quat(mp);
    flatness_map = A01_FlatnessMap;
    fprintf('Selected model: A01 (13-State Quaternion representation)\n');
else
    quad = A00_SimplifiedModel(mp);
    flatness_map = A00_FlatnessMap;
    fprintf('Selected model: A00 (12-State Euler Z-Y-X representation)\n');
end

fprintf('OpenLoopLemniscate \n');
fprintf('Drone: m=%.3f kg, T_hover=%.3f N\n\n', mp.m, mp.m*mp.g);

%% TRAJECTORY CONFIGURATION
% Select Lemniscate formulation (Gerono)
lp = tp.Lemniscate2;

traj = struct( ...
    'name', 'Lemniscate2', ...
    'obj',  ShapeLemniscate2(lp), ...
    'col',  [0.6 0.2 0.8]); 

T_final = 1 * tp.T_duration_lemniscate2;
fprintf('Simulation time: %.3f s\n\n', T_final);

%% ODE OPTIONS
opts_ode = odeset('RelTol', 1e-7, 'AbsTol', 1e-8);

%% INITIAL CONDITION FROM FLATNESS MAP
[s0, ds0, dds0, ddds0, dddds0] = traj.obj.get_flat_outputs(0);
[x0, u0] = flatness_map.map(mp, s0, ds0, dds0, ddds0, dddds0);

% Extract Euler equivalent for standard command-line logging
if use_quaternions
    [init_phi, init_theta, init_psi] = Quat_to_euler(x0(7:10));
    init_p = x0(11); init_q = x0(12); init_r = x0(13);
else
    init_phi = x0(7); init_theta = x0(8); init_psi = x0(9);
    init_p = x0(10); init_q = x0(11); init_r = x0(12);
end

fprintf('Initial condition:\n');
fprintf('  x_I   = %.2f m\n', x0(1));
fprintf('  y_I   = %.2f m\n', x0(2));
fprintf('  z_I   = %.2f m\n', x0(3));
fprintf('  vx_B  = %.2f m/s\n', x0(4));
fprintf('  vy_B  = %.2f m/s\n', x0(5));
fprintf('  vz_B  = %.2f m/s\n', x0(6));
fprintf('  phi   = %.2f deg\n', init_phi*180/pi);
fprintf('  theta = %.2f deg\n', init_theta*180/pi);
fprintf('  psi   = %.2f deg\n', init_psi*180/pi);
fprintf('  p     = %.2f rad/s\n', init_p);
fprintf('  q     = %.2f rad/s\n', init_q);
fprintf('  r     = %.2f rad/s\n\n', init_r);

fprintf('Initial inputs:\n');
fprintf('  T = %.4f N\n', u0(1));
fprintf('  L = %.3e Nm\n', u0(2));
fprintf('  M = %.3e Nm\n', u0(3));
fprintf('  N = %.3e Nm\n\n', u0(4));

%% ODE INTEGRATION
[t_sol, x_sol] = ode45( ...
    @(t, x) ode_dynamics(t, x, quad, traj.obj, mp, use_quaternions, flatness_map), ...
    [0 T_final], x0, opts_ode);

fprintf('Final integration time: %.4f s\n\n', t_sol(end));

%% RECONSTRUCT REFERENCES
N = length(t_sol);
if use_quaternions
    x_ref = zeros(N, 13);
else
    x_ref = zeros(N, 12);
end
u_ref = zeros(N, 4);

for k = 1:N
    [s, ds, dds, ddds, dddds] = traj.obj.get_flat_outputs(t_sol(k));
    [xr, ur] = flatness_map.map(mp, s, ds, dds, ddds, dddds);
    x_ref(k, :) = xr';
    u_ref(k, :) = ur';
end

% Convert simulated and reference trajectories to 12-state (Euler) for plotting 
x_sol_euler = zeros(N, 12);
x_ref_euler = zeros(N, 12);

x_sol_euler(:, 1:6) = x_sol(:, 1:6);
x_ref_euler(:, 1:6) = x_ref(:, 1:6);

if use_quaternions
    x_sol_euler(:, 10:12) = x_sol(:, 11:13);
    x_ref_euler(:, 10:12) = x_ref(:, 11:13);
    for k = 1:N
        [phi_s, theta_s, psi_s] = Quat_to_euler(x_sol(k, 7:10)');
        x_sol_euler(k, 7:9) = [phi_s, theta_s, psi_s];
        
        [phi_r, theta_r, psi_r] = Quat_to_euler(x_ref(k, 7:10)');
        x_ref_euler(k, 7:9) = [phi_r, theta_r, psi_r];
    end
else
    x_sol_euler(:, 7:12) = x_sol(:, 7:12);
    x_ref_euler(:, 7:12) = x_ref(:, 7:12);
end

% Unwrap yaw angle to handle rotation discontinuities
x_ref_euler(:, 9) = unwrap(x_ref_euler(:, 9));
x_sol_euler(:, 9) = unwrap(x_sol_euler(:, 9));

%% POSITION ERROR
pos_err = vecnorm(x_sol_euler(:, 1:3) - x_ref_euler(:, 1:3), 2, 2);
fprintf('Max position error: %.3e m\n', max(pos_err));

%% FIGURE 1 — STATE TRACKING
state_names = {
    'X_I [m]','Y_I [m]','Z_{NED} [m]',...
    'u_b [m/s]','v_b [m/s]','w_b [m/s]',...
    '\phi [deg]','\theta [deg]','\psi [deg]',...
    'p [rad/s]','q [rad/s]','r [rad/s]'};

angle_scale = [1 1 1 1 1 1 180/pi 180/pi 180/pi 1 1 1];

figure('Name', 'Lemniscate State Tracking', 'Color', 'w', 'Position', [30 50 1600 850]);
tiledlayout(4, 3, 'Padding', 'compact', 'TileSpacing', 'compact');

for j = 1:12
    nexttile;
    plot(t_sol, x_ref_euler(:, j)*angle_scale(j), 'r--', 'LineWidth', 1.8); hold on;
    plot(t_sol, x_sol_euler(:, j)*angle_scale(j), 'b', 'LineWidth', 1.2);
    ylabel(state_names{j}); grid on;

    if j <= 3, title(state_names{j}); end
    if j > 9,  xlabel('t [s]'); end
    if j == 1, legend('ref', 'sim', 'Location', 'best'); end
end
sgtitle(sprintf('Lemniscate tracking  |  max error = %.2e m', max(pos_err)));

%% FIGURE 2 — INPUTS
figure('Name', 'Lemniscate Inputs', 'Color', 'w', 'Position', [30 50 900 500]);
input_labels = {'T [N]', 'L [N m]', 'M [N m]', 'N [N m]'};
T_hover = mp.m * mp.g;

for j = 1:4
    subplot(4, 1, j);
    plot(t_sol, u_ref(:, j), 'Color', traj.col, 'LineWidth', 1.8); grid on;
    if j == 1
        yline(T_hover, 'k--', 'T_{hover}');
    else
        yline(0, 'k--');
    end
    ylabel(input_labels{j});
end
xlabel('t [s]'); sgtitle('Reference inputs');

%% FIGURE 3 — POSITION ERROR
figure('Name', 'Position Error', 'Color', 'w');
plot(t_sol, pos_err, 'LineWidth', 2); grid on;
xlabel('t [s]'); ylabel('||e_p|| [m]');
title('Position tracking error');

%% FIGURE 4 — 3D DRONE ANIMATION 
fig3 = figure('Name', '3D Lemniscate Animation', 'Color', 'w', 'Position', [900 100 800 750]);
ax3 = axes(fig3);
PlotDrone.setup_axes(ax3, '3D Drone — Lemniscate Trajectory');
hold(ax3, 'on');

sz = 3.0; 
PlotDrone.draw_ground(ax3, [-sz sz], [-sz sz], 0);

plot3(ax3, x_ref_euler(:, 1), x_ref_euler(:, 2), x_ref_euler(:, 3), 'k--', 'LineWidth', 1.5);
h_traj = plot3(ax3, NaN, NaN, NaN, 'Color', traj.col, 'LineWidth', 1.5);
h_drone = PlotDrone.init(ax3, mp);

pos0 = x_sol_euler(1, 1:3)';
R0 = Euler2R_ZYX(x_sol_euler(1, 7), x_sol_euler(1, 8), x_sol_euler(1, 9));
PlotDrone.update(h_drone, pos0, R0);

x_all = [x_ref_euler(:, 1); x_sol_euler(:, 1)];
y_all = [x_ref_euler(:, 2); x_sol_euler(:, 2)];
z_all = [x_ref_euler(:, 3); x_sol_euler(:, 3)];
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
    pos_k = x_sol_euler(k, 1:3)';
    phi_k = x_sol_euler(k, 7);
    the_k = x_sol_euler(k, 8);
    psi_k = x_sol_euler(k, 9);
    R_k   = Euler2R_ZYX(phi_k, the_k, psi_k);

    PlotDrone.update(h_drone, pos_k, R_k);
    set(h_traj, 'XData', x_sol_euler(1:k, 1), 'YData', x_sol_euler(1:k, 2), 'ZData', x_sol_euler(1:k, 3));

    title(ax3, sprintf([ ...
        'Lemniscate  |  t = %.2f s  |  ', ...
        '\\phi = %.1f°  ', ...
        '\\theta = %.1f°  ', ...
        '\\psi = %.1f°'], ...
        t_sol(k), phi_k*180/pi, the_k*180/pi, psi_k*180/pi));
    drawnow;
end

fprintf('\n End simulation \n');

%% LOCAL FUNCTIONS
function dxdt = ode_dynamics(t, x, quad, traj_obj, mp, use_quaternions, flatness_map)
    % Reconstruct trajectory state and derivatives
    [s, ds, dds, ddds, dddds] = traj_obj.get_flat_outputs(t);

    % Get control input commands from Flatness Map
    [~, uref] = flatness_map.map(mp, s, ds, dds, ddds, dddds);

    % Renormalize quaternion state inside ODE loop to avoid integration drift 
    if use_quaternions
        q = x(7:10);
        x(7:10) = q / norm(q);
    end

    % Evaluate dynamics
    dxdt = quad.dynamics(t, x, uref);
end
