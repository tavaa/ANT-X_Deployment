%% Test open-loop for the quadcopter Degrees of Freedom (DoF) -- B01 model.
%
% The 4 tests isolate a single DOF at a time, verifying that the flatness
% map produces physically correct motor thrust inputs and that the
% nonlinear model tracks the reference.
%
% Works with quad.delay_motors mode:
%   true  : 16-state plant with motor delay.
%   false : 12-state, no-delay plant.
%
% B01_FlatnessMap.map ALWAYS returns a 16-element state (rigid body +
% no-delay T1..T4 padding) regardless of quad.delay_motors -- this
% script truncates to n_states = 12+4*quad.delay_motors before handing
% anything to ode45, to match whichever model quad actually implements.
%
%  TEST  1 -- PITCH
%       Expected:  theta varies, phi=0, psi=0
%       FRD:       theta > 0 = nose up
%
%  TEST  2 -- ROLL
%       Expected:  phi varies, theta=0, psi=0
%       FRD:       phi > 0 = right wing down
%
%  TEST  3 -- YAW
%       Expected:  psi varies, phi=theta=0,  T ~= mg
%
%  TEST  4 -- VERTICAL (altitude)
%       Expected:  z varies, phi=theta=psi=0,  T1=T2=T3=T4 (pure collective)
%
%  Conventions:
%   Inertial frame : NED  (x=North, y=East, z=Down positive)
%   Body frame     : FRD  (x=Forward, y=Right, z=Down)
%   Velocity state : NED (v_n, v_e, v_d) -- B01 is natively inertial-velocity.

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

quad = B01_ANTX_quadcopter(mp, true);    % delay_motors=false: 12-state, no-delay plant
flatness_map = B01_FlatnessMap;

n_states = 12 + 4*quad.delay_motors;      % 16 or 12, matches quad exactly

fprintf('Selected model: B01 (delay_motors=%d, %d states)\n', quad.delay_motors, n_states);
fprintf('Drone: m=%.3f kg,  T_hover_total=%.3f N,  T_hover_per_motor=%.4f N\n\n', ...
    mp.m, mp.m*mp.g, mp.m*mp.g/4);

%% ODE OPTIONS
opts_ode = odeset('RelTol', 1e-7, 'AbsTol', 1e-8);

%% TEST DEFINITIONS
p1 = tp.Pitch;   t1 = ShapePitchTest(p1);   T1 = tp.T_duration_pitch;
p2 = tp.Roll;    t2 = ShapeRollTest(p2);    T2 = tp.T_duration_roll;
p3 = tp.Yaw;     t3 = ShapeYawTest(p3);     T3 = tp.T_duration_yaw;
p4 = tp.Thrust;  t4 = ShapeThrustTest(p4);  T4 = tp.T_duration_thrust;

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

    % Extract initial condition from Flatness Map, truncate to n_states
    [s0, ds0, dds0, ddds0, dddds0] = tr.obj.get_flat_outputs(0);
    [x0_full, u0] = flatness_map.map(mp, s0, ds0, dds0, ddds0, dddds0);
    x0 = x0_full(1:n_states);

    init_phi = x0(7); init_theta = x0(8); init_psi = x0(9);
    init_p   = x0(10); init_q = x0(11); init_r = x0(12);

    fprintf('  Initial Position (NED): x_I=%.2f  y_I=%.2f  z_I=%.2f\n', x0(1), x0(2), x0(3));
    fprintf('  Initial Inertial Velocities (NED): vn=%.2f  ve=%.2f  vd=%.2f\n', x0(4), x0(5), x0(6));
    fprintf('  Initial Euler Angles: phi=%.2f  theta=%.2f  psi=%.2f deg\n', ...
        init_phi*180/pi, init_theta*180/pi, init_psi*180/pi);
    fprintf('  Initial Body Rates: p=%.3f  q=%.3f  r=%.3f rad/s\n', init_p, init_q, init_r);
    if quad.delay_motors
        fprintf('  Initial Motor Thrusts (x0): T1=%.4f  T2=%.4f  T3=%.4f  T4=%.4f\n', ...
            x0(13), x0(14), x0(15), x0(16));
    end
    fprintf('  Initial Inputs (u0): T1=%.4f  T2=%.4f  T3=%.4f  T4=%.4f\n', ...
        u0(1), u0(2), u0(3), u0(4));

    % ODE45 Integration
    try
        [t_sol, x_sol] = ode45(...
            @(t, x) ode_dynamics(t, x, quad, tr.obj, mp, flatness_map), ...
            [0, tr.T], x0, opts_ode);
    catch ME
        fprintf('  ERROR: %s\n\n', ME.message);
        results{i} = struct('name', tr.name, 'failed', true, 'col', tr.col);
        continue;
    end

    % Reference sampling from Flatness Map, truncated to n_states
    N_steps = length(t_sol);
    x_ref = zeros(N_steps, n_states);
    u_ref = zeros(N_steps, 4);

    for k = 1:N_steps
        [s, ds, dds, ddds, dddds] = tr.obj.get_flat_outputs(t_sol(k));
        [xr_full, ur] = flatness_map.map(mp, s, ds, dds, ddds, dddds);
        x_ref(k, :) = xr_full(1:n_states)';
        u_ref(k, :) = ur';
    end

    % Position tracking error
    pos_err = vecnorm(x_sol(:, 1:3) - x_ref(:, 1:3), 2, 2);

    % Motor thrust tracking error -- only meaningful if delay_motors=true
    if quad.delay_motors
        thrust_err = vecnorm(x_sol(:, 13:16) - x_ref(:, 13:16), 2, 2);
    else
        thrust_err = [];
    end

    fprintf('  Steps: %d,  max pos error: %.2e m', N_steps, max(pos_err));
    if quad.delay_motors
        fprintf(',  max thrust-state error: %.2e N', max(thrust_err));
    end
    fprintf('\n');
    fprintf('  max |phi|=%.2f  |theta|=%.2f  |psi|=%.2f deg\n', ...
        max(abs(x_ref(:, 7)))*180/pi, max(abs(x_ref(:, 8)))*180/pi, max(abs(x_ref(:, 9)))*180/pi);
    fprintf('  max |T1|=%.4f  |T2|=%.4f  |T3|=%.4f  |T4|=%.4f\n\n', ...
        max(abs(u_ref(:, 1))), max(abs(u_ref(:, 2))), max(abs(u_ref(:, 3))), max(abs(u_ref(:, 4))));

    results{i} = struct('name', tr.name, 'failed', false, 'col', tr.col, ...
        't', t_sol, 'x', x_sol, 'x_ref', x_ref, 'u_ref', u_ref, ...
        'err', pos_err, 'thrust_err', thrust_err, 'params', tr.obj.params);
end

%% GRAPHS, PLOTS AND INSIGHTS
state_names = {
    'r_n [m]','r_e [m]','r_d [m]',...
    'v_n [m/s]','v_e [m/s]','v_d [m/s]',...
    '\phi [deg]','\theta [deg]','\psi [deg]',...
    'p [rad/s]','q [rad/s]','r [rad/s]'};

angle_scale = [1 1 1 1 1 1 180/pi 180/pi 180/pi 1 1 1];
T_hover_pm = mp.m * mp.g / 4;

for i = 1:n_tests
    r = results{i};
    if r.failed, continue; end

    %% FIG 1: Rigid-body State Tracking (states 1-12, always present)
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

    %% FIG 1b: Motor Thrust States (13-16) -- only if delay_motors=true
    if quad.delay_motors
        fig1b = figure('Name', ['Motor thrusts: ' r.name], 'Color', 'w', 'Position', [30 570 900 350]);
        motor_labels = {'T_1 [N]', 'T_2 [N]', 'T_3 [N]', 'T_4 [N]'};
        for j = 1:4
            subplot(1, 4, j);
            plot(r.t, r.x_ref(:, 12+j), 'r--', 'LineWidth', 1.8); hold on;
            plot(r.t, r.x(:, 12+j),     'b',   'LineWidth', 1.2);
            yline(T_hover_pm, 'k:');
            ylabel(motor_labels{j}); xlabel('t [s]'); grid on;
            if j == 1, legend('ref (no delay)', 'sim (with delay)', 'Location', 'best'); end
        end
        sgtitle(['Motor thrust tracking: ' r.name ...
            sprintf('  (max err = %.2e N)', max(r.thrust_err))]);
    end

    %% FIG 2: Control Inputs (desired motor thrusts)
    fig2 = figure('Name', ['Inputs: ' r.name], 'Color', 'w', 'Position', [30 950 900 500]);
    input_labels = {'T_1^d [N]', 'T_2^d [N]', 'T_3^d [N]', 'T_4^d [N]'};
    for j = 1:4
        subplot(4, 1, j);
        plot(r.t, r.u_ref(:, j), 'Color', r.col, 'LineWidth', 1.8);
        hold on;
        y_data = r.u_ref(:, j);
        y_min = min(y_data); y_max = max(y_data);
        y_span = max(y_max - y_min, 1e-6);
        y_margin = 0.15 * y_span;
        ylim([y_min - y_margin, max(y_max, T_hover_pm) + y_margin]);
        xlim([r.t(1), r.t(end)]);
        yline(T_hover_pm, 'k--', 'T_{hover}/4', 'LabelHorizontalAlignment', 'left');
        ylabel(input_labels{j}); grid on;
    end
    xlabel('t [s]');
    sgtitle(['Reference inputs (desired motor thrusts): ' r.name]);

    %% FIG 3: 3D Animation
    fig3 = figure('Name', ['Anim: ' r.name], 'Color', 'w', 'Position', [900 100 800 750]);
    ax3  = axes(fig3);
    PlotDrone2.setup_axes(ax3, ['3D Drone -- ' r.name]);
    hold(ax3, 'on');

    sz = 1.5;
    PlotDrone2.draw_ground(ax3, [-sz sz], [-sz sz], 0);

    plot3(ax3, r.x_ref(:, 1), r.x_ref(:, 2), r.x_ref(:, 3), 'k--', 'LineWidth', 1.5);
    h_traj = plot3(ax3, NaN, NaN, NaN, 'Color', r.col, 'LineWidth', 1.2);
    h_drone = PlotDrone2.init(ax3, mp);

    pos0 = r.x(1, 1:3)';
    R0   = Euler2R_ZYX(r.x(1, 7), r.x(1, 8), r.x(1, 9));
    PlotDrone2.update(h_drone, pos0, R0);

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

        PlotDrone2.update(h_drone, pos_k, R_k);
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

    sgtitle(['Euler angles -- ' r.name]);
end

%% TEST COMPARISON -- desired motor thrusts across all tests
figure('Name', 'Input summary', 'Color', 'w', 'Position', [50 50 1200 500]);
input_labels2 = {'T_1^d [N]', 'T_2^d [N]', 'T_3^d [N]', 'T_4^d [N]'};
for j = 1:4
    subplot(1, 4, j);
    hold on; grid on;
    for i = 1:n_tests
        r = results{i};
        if r.failed, continue; end
        plot(r.t, r.u_ref(:, j), 'Color', r.col, 'LineWidth', 1.5, 'DisplayName', r.name);
    end
    yline(T_hover_pm, 'k--', 'T_{hover}/4');
    if j == 1
        legend('Location', 'best', 'FontSize', 7);
    end
    title(['Input ' input_labels2{j}]);
    xlabel('t [s]'); ylabel(input_labels2{j});
end
sgtitle('Reference inputs (motor thrusts), all tests');

fprintf('End Simulation \n');

%% LOCAL FUNCTIONS
function dxdt = ode_dynamics(t, x, quad, traj_obj, mp, flatness_map)
    % ODE_DYNAMICS  RHS for ode45: sample the reference trajectory at
    % time t, reconstruct the desired motor thrusts via flatness, and
    % evaluate the plant dynamics at the current state x under that
    % open-loop input.
    [s, ds, dds, ddds, dddds] = traj_obj.get_flat_outputs(t);

    % uref is always 4-elements and means the same thing (T_i^d) in both
    % delay_motors modes -- no truncation needed on the input side.
    [~, uref] = flatness_map.map(mp, s, ds, dds, ddds, dddds);

    % Evaluate dynamics (12 or 16 states, matching quad.delay_motors)
    dxdt = quad.dynamics(t, x, uref);
end