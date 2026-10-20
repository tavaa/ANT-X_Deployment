%% Test open-loop for the quadcopter Degrees of Freedom (DoF).
%
% The 4 tests isolate a single DOF at a time, verifying that the flatness
% map produces physically correct inputs and that the nonlinear model
% tracks the reference with machine-level error (for small angles).
%
%  TEST  1 — PITCH
%       Expected:  theta varies, phi=0, psi=0
%                  M ≠ 0,  L=N=0,  T varies
%       FRD:       theta > 0 = nose up
%
%  TEST  2 — ROLL
%       Expected:  phi varies, theta=0, psi=0
%                  L ≠ 0,  M=N=0,  T varies
%       FRD:       phi > 0 = right wing down 
%
%  TEST  3 — YAW
%       Expected:  psi varies, phi=theta=0
%                  N ≠ 0,  L=M=0,  T = mg (constant)
%       Note:      N ≠ 0 because psi_ddot ≠ 0 (angular acceleration)
%
%  TEST  4 — VERTICAL (altitude)
%       Expected:  z varies, phi=theta=psi=0
%                  L=M=N=0,  T varies around mg
%
%  Conventions:
%   Inertial frame : NED  (x=North, y=East, z=Down positive)
%   Body frame     : FRD  (x=Forward, y=Right, z=Down)

clear; clc; close all;

%% PATHS
addpath(fullfile('models'));
addpath(fullfile('flatness'));
addpath(fullfile('configs'));
addpath(fullfile('trajectories'));
addpath(fullfile('utils'));

%% SETTINGS
use_quaternions = false; % Set to true for A01 (quaternion), false for A00 (Euler Z-Y-X)

%% MODEL SELECTION
mp = ModelParameters;
tp = TrajectoryParameters;

if use_quaternions
    quad = A01_SimplifiedModel_quat(mp);
    flatness_map = A01_FlatnessMap;
    fprintf('Selected model: A01 (13-State Quaternion representation)\n');
else
    %quad = A00_SimplifiedModel(mp);
    %flatness_map = A00_FlatnessMap;
    quad = A00_SimplifiedModel(mp);
    flatness_map = A00_FlatnessMap;
    fprintf('Selected model: A00 (12-State Euler Z-Y-X representation)\n');
end

fprintf('Drone: m=%.3f kg,  T_hover=%.3f N\n\n', mp.m, mp.m*mp.g);

%% ODE OPTIONS
opts_ode = odeset('RelTol', 1e-7, 'AbsTol', 1e-8);

%% TEST DEFINITIONS
% PITCH 
p1 = tp.Pitch;
t1 = ShapePitchTest(p1);
T1 = tp.T_duration_pitch;

% ROLL 
p2 = tp.Roll;
t2 = ShapeRollTest(p2);
T2 = tp.T_duration_roll;

% YAW 
p3 = tp.Yaw;
t3 = ShapeYawTest(p3);
T3 = tp.T_duration_yaw;

% VERTICAL
p4 = tp.Thrust;
t4 = ShapeThrustTest(p4);
T4 = tp.T_duration_thrust;

% Collect tests
tests = {
    struct('name','PitchTest','obj',t1,'T',T1,'col',[0.2 0.5 0.8])
    struct('name','RollTest', 'obj',t2,'T',T2,'col',[0.8 0.2 0.2])
    struct('name','YawTest',  'obj',t3,'T',T3,'col',[0.1 0.7 0.3])
    struct('name','VertTest', 'obj',t4,'T',T4,'col',[0.7 0.4 0.0])
};

n_tests = numel(tests);
results = cell(n_tests, 1);

%% SIMULATION LOOP
for i = 1:n_tests
    tr = tests{i};
    fprintf('%s: \n', tr.name);

    % Extract initial condition from Flatness Map
    [s0, ds0, dds0, ddds0, dddds0] = tr.obj.get_flat_outputs(0);
    [x0, u0] = flatness_map.map(mp, s0, ds0, dds0, ddds0, dddds0);
    
    % Reconstruct initial Euler angles for clean logging
    if use_quaternions
        [init_phi, init_theta, init_psi] = Quat_to_euler(x0(7:10));
        init_p = x0(11); init_q = x0(12); init_r = x0(13);
    else
        init_phi = x0(7); init_theta = x0(8); init_psi = x0(9);
        init_p = x0(10); init_q = x0(11); init_r = x0(12);
    end

    fprintf('  Initial Position (NED): x_I=%.2f  y_I=%.2f  z_I=%.2f\n', x0(1), x0(2), x0(3));
    fprintf('  Initial Body Velocities (FRD): vx_B=%.2f  vy_B=%.2f  vz_B=%.2f\n', x0(4), x0(5), x0(6));
    fprintf('  Initial Euler Angles: phi=%.2f  theta=%.2f  psi=%.2f deg\n', ...
        init_phi*180/pi, init_theta*180/pi, init_psi*180/pi);
    fprintf('  Initial Body Rates: p=%.3f  q=%.3f  r=%.3f rad/s\n', init_p, init_q, init_r);
    fprintf('  Initial Inputs (u0): T=%.4f  L=%.2e  M=%.2e  N=%.2e\n', u0(1), u0(2), u0(3), u0(4));

    % ODE45 Integration
    try
        [t_sol, x_sol] = ode45(...
            @(t, x) ode_dynamics(t, x, quad, tr.obj, mp, use_quaternions, flatness_map), ...
            [0, tr.T], x0, opts_ode);
    catch ME
        fprintf('  ERROR: %s\n\n', ME.message);
        results{i} = struct('name', tr.name, 'failed', true, 'col', tr.col);
        continue;
    end

    % Reference sampling from Flatness Map
    N_steps = length(t_sol);
    if use_quaternions
        x_ref = zeros(N_steps, 13);
    else
        x_ref = zeros(N_steps, 12);
    end
    u_ref = zeros(N_steps, 4);
    
    for k = 1:N_steps
        [s, ds, dds, ddds, dddds] = tr.obj.get_flat_outputs(t_sol(k));
        [xr, ur] = flatness_map.map(mp, s, ds, dds, ddds, dddds);
        x_ref(k, :) = xr';
        u_ref(k, :) = ur';
    end
    
    % Position tracking error check
    pos_err = vecnorm(x_sol(:, 1:3) - x_ref(:, 1:3), 2, 2);

    % Dynamic state conversion to Euler-equivalent representation for plots and visual
    x_sol_euler = zeros(N_steps, 12);
    x_ref_euler = zeros(N_steps, 12);
    
    x_sol_euler(:, 1:6) = x_sol(:, 1:6);
    x_ref_euler(:, 1:6) = x_ref(:, 1:6);
    
    if use_quaternions
        x_sol_euler(:, 10:12) = x_sol(:, 11:13);
        x_ref_euler(:, 10:12) = x_ref(:, 11:13);
        for k = 1:N_steps
            [phi_s, theta_s, psi_s] = Quat_to_euler(x_sol(k, 7:10)');
            x_sol_euler(k, 7:9) = [phi_s, theta_s, psi_s];
            
            [phi_r, theta_r, psi_r] = Quat_to_euler(x_ref(k, 7:10)');
            x_ref_euler(k, 7:9) = [phi_r, theta_r, psi_r];
        end
    else
        x_sol_euler(:, 7:12) = x_sol(:, 7:12);
        x_ref_euler(:, 7:12) = x_ref(:, 7:12);
    end

    fprintf('  Steps: %d,  max pos error: %.2e m\n', N_steps, max(pos_err));
    fprintf('  max |phi|=%.2f  |theta|=%.2f  |psi|=%.2f deg\n', ...
        max(abs(x_ref_euler(:, 7)))*180/pi, max(abs(x_ref_euler(:, 8)))*180/pi, max(abs(x_ref_euler(:, 9)))*180/pi);
    fprintf('  max |T|=%.4f  |L|=%.2e  |M|=%.2e  |N|=%.2e\n\n', ...
        max(abs(u_ref(:, 1))), max(abs(u_ref(:, 2))), max(abs(u_ref(:, 3))), max(abs(u_ref(:, 4))));

    results{i} = struct('name', tr.name, 'failed', false, 'col', tr.col,...
        't', t_sol, 'x', x_sol_euler, 'x_ref', x_ref_euler, 'u_ref', u_ref, 'err', pos_err,...
        'params', tr.obj.params);
end

%% GRAPHS, PLOTS AND INSIGHTS
state_names = {
    'X_I [m]','Y_I [m]','Z_{NED} [m]',...
    'xb [m/s]','yb [m/s]','zb [m/s]',...
    '\phi [deg]','\theta [deg]','\psi [deg]',...
    'p [rad/s]','q [rad/s]','r [rad/s]'};

angle_scale = [1 1 1 1 1 1 180/pi 180/pi 180/pi 1 1 1]; 

for i = 1:n_tests
    r = results{i};
    if r.failed, continue; end

    %% FIG 1: State Tracking
    fig1 = figure('Name', ['State: ' r.name], 'Color', 'w', 'Position', [30 50 900 500]);
    tiledlayout(4, 3, 'Padding', 'compact', 'TileSpacing', 'compact');

    for j = 1:12
        nexttile;
        plot(r.t, r.x_ref(:, j)*angle_scale(j), 'r--', 'LineWidth', 1.8); hold on;
        plot(r.t, r.x(:, j)*angle_scale(j),     'b',   'LineWidth', 1.2);
        ylabel(state_names{j}); grid on;
        if j <= 3, title(state_names{j}); end
        if j > 9, xlabel('t [s]'); end
        if j == 1, legend('ref', 'sim', 'Location', 'best'); end
    end
    sgtitle(['State tracking: ' r.name sprintf('  (max pos err = %.2e m)', max(r.err))]);

    %% FIG 2: Control Inputs
    fig2 = figure('Name', ['Inputs: ' r.name], 'Color', 'w', 'Position', [30 50 900 500]);
    input_labels = {'T [N]', 'L [N\cdotm]', 'M [N\cdotm]', 'N [N\cdotm]'};
    T_hover = mp.m * mp.g;
    for j = 1:4
        subplot(4, 1, j);
        plot(r.t, r.u_ref(:, j), 'Color', r.col, 'LineWidth', 1.8);
        if j == 1, yline(T_hover, 'k--', 'T_{hover}', 'LabelHorizontalAlignment', 'right'); end
        if j > 1,  yline(0, 'k--'); end
        ylabel(input_labels{j}); grid on;
    end
    xlabel('t [s]');
    sgtitle(['Reference inputs: ' r.name]);

    %% FIG 3: 3D Animation [1]
    fig3 = figure('Name', ['Anim: ' r.name], 'Color', 'w', 'Position', [900 100 800 750]);
    ax3  = axes(fig3);
    PlotDrone.setup_axes(ax3, ['3D Drone — ' r.name]);
    hold(ax3, 'on');

    sz = 1.5;
    PlotDrone.draw_ground(ax3, [-sz sz], [-sz sz], 0);

    plot3(ax3, r.x_ref(:, 1), r.x_ref(:, 2), r.x_ref(:, 3), 'k--', 'LineWidth', 1.5);
    h_traj = plot3(ax3, NaN, NaN, NaN, 'Color', r.col, 'LineWidth', 1.2);
    h_drone = PlotDrone.init(ax3, mp);

    pos0 = r.x(1, 1:3)';
    R0   = Euler2R_ZYX(r.x(1, 7), r.x(1, 8), r.x(1, 9));
    PlotDrone.update(h_drone, pos0, R0);

    x_all = [r.x_ref(:, 1); r.x(:, 1)];
    y_all = [r.x_ref(:, 2); r.x(:, 2)];
    z_all = [r.x_ref(:, 3); r.x(:, 3)];
    margin = 0.3;
    xlim(ax3, [min(x_all)-margin, max(x_all)+margin]);
    ylim(ax3, [min(y_all)-margin, max(y_all)+margin]);
    zlim(ax3, [min(z_all)-margin, max(z_all)+margin]);

    annotation(fig3, 'textbox', [0.01 0.01 0.28 0.14],...
        'String', {'FRD body frame:', '  x_B (red)   = Forward', '  y_B (green) = Right', '  z_B (blue)  = Down',...
                   '  \phi > 0 : right wing DOWN', '  \theta > 0 : nose UP', '  \psi > 0 : CW rotation'},...
        'FontSize', 7, 'BackgroundColor', 'w', 'EdgeColor', [0.5 0.5 0.5]);

    skip = max(1, floor(length(r.t)/120));
    for k = 1:skip:length(r.t)
        pos_k = r.x(k, 1:3)';
        phi_k = r.x(k, 7);
        the_k = r.x(k, 8);
        psi_k = r.x(k, 9);
        R_k   = Euler2R_ZYX(phi_k, the_k, psi_k);

        PlotDrone.update(h_drone, pos_k, R_k);
        set(h_traj, 'XData', r.x(1:k, 1), 'YData', r.x(1:k, 2), 'ZData', r.x(1:k, 3));
        title(ax3, sprintf('%s  |  t=%.2f s  |  \\phi=%.1f°  \\theta=%.1f°  \\psi=%.1f°', ...
            r.name, r.t(k), phi_k*180/pi, the_k*180/pi, psi_k*180/pi));
        drawnow;
    end

    %% FIG 4: Euler Angles
    fig4 = figure('Name', ['Angles: ' r.name], 'Color', 'w', 'Position', [30 50 900 350]);
    subplot(1, 3, 1);
    plot(r.t, r.x_ref(:, 7)*180/pi, 'r--', 'LineWidth', 1.8); hold on;
    plot(r.t, r.x(:, 7)*180/pi, 'b', 'LineWidth', 1.2);
    yline(0, 'k:'); xlabel('t [s]'); ylabel('\phi [deg]');
    title('\phi  (roll)'); grid on; legend('ref', 'sim');

    subplot(1, 3, 2);
    plot(r.t, r.x_ref(:, 8)*180/pi, 'r--', 'LineWidth', 1.8); hold on;
    plot(r.t, r.x(:, 8)*180/pi, 'b', 'LineWidth', 1.2);
    yline(0, 'k:'); xlabel('t [s]'); ylabel('\theta [deg]');
    title('\theta  (pitch)'); grid on;

    subplot(1, 3, 3);
    plot(r.t, r.x_ref(:, 9)*180/pi, 'r--', 'LineWidth', 1.8); hold on;
    plot(r.t, r.x(:, 9)*180/pi, 'b', 'LineWidth', 1.2);
    yline(0, 'k:'); xlabel('t [s]'); ylabel('\psi [deg]');
    title('\psi  (yaw)'); grid on;

    sgtitle(['Euler angles — ' r.name]);
end

%% TEST COMPARISON
figure('Name', 'Input summary', 'Color', 'w', 'Position', [50 50 1200 500]);
input_labels2 = {'T [N]', 'L [N\cdotm]', 'M [N\cdotm]', 'N [N\cdotm]'};
for j = 1:4
    subplot(1, 4, j);
    hold on; grid on;
    for i = 1:n_tests
        r = results{i};
        if r.failed, continue; end
        plot(r.t, r.u_ref(:, j), 'Color', r.col, 'LineWidth', 1.5, 'DisplayName', r.name);
    end
    if j == 1
        yline(T_hover, 'k--', 'T_{hover}');
        legend('Location', 'best', 'FontSize', 7);
        title('T (hover)');
    else
        yline(0, 'k--');
        title(['Input ' input_labels2{j}]);
    end
    xlabel('t [s]'); ylabel(input_labels2{j});
end
sgtitle('Reference inputs all tests');

fprintf('End Simulation \n');

%% LOCAL FUNCTIONS
function dxdt = ode_dynamics(t, x, quad, traj_obj, mp, use_quaternions, flatness_map)
    % Sample trajectory state and derivatives
    [s, ds, dds, ddds, dddds] = traj_obj.get_flat_outputs(t);

    % Get reference inputs from Flatness Map
    [~, uref] = flatness_map.map(mp, s, ds, dds, ddds, dddds);

    % Renormalize quaternion state inside ODE loop to avoid drift
    if use_quaternions
        q = x(7:10);
        x(7:10) = q / norm(q);
    end

    % Evaluate dynamics
    dxdt = quad.dynamics(t, x, uref);
end

