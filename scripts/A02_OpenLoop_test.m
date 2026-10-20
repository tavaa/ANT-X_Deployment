%% Test open-loop for the quadcopter — A02 model,
%
%  Circle / Lemniscate / Spiral / Lissajous trajectory tracking using differential
%  flatness. Same test repeated for all the shapes in a single script.
%
%  MATLAB only OpenLoop test:
%  for deployment use deployment/OpenLoop_Test_Shapes.m
%
%  Frames:
%       Inertial : NED
%       Body     : FRD
%       Velocity state : NED (v_x, v_y, v_z).
%
%  Model: A02_SimplifiedModel_Deployment.

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

quad = A02_SimplifiedModel_Deployment(mp);
flatness_map = A02_FlatnessMap;
fprintf('Selected model: A02 (12-states model)\n');

fprintf('OpenLoop A02 — Circle / Lemniscate / Spiral / Lissajous\n');
fprintf('Drone: m=%.3f kg, T_hover_total=%.3f N', ...
    mp.m, mp.m*mp.g);

%% TRAJECTORY LIST
% OPTIONAL: test also Circle Shape Deployment
trajs = { ...
    struct('name', 'Circle',      'obj', ShapeCircle(tp.Circle),           'col', [0.1 0.45 0.9], 'T_final', 1*tp.T_duration_circle,       'ground_size', 2.5), ...
    struct('name', 'Circle Deployment',      'obj', ShapeCircleDeployment(tp.CircleDeployment), 'col', [0.1 0.45 0.9], 'T_final', 1*tp.T_duration_circle_deployment, 'ground_size', 2.5), ...
    struct('name', 'Lemniscate', 'obj', ShapeLemniscate(tp.Lemniscate), 'col', [0.6 0.2  0.8], 'T_final', 1*tp.T_duration_lemniscate,  'ground_size', 3.0), ...
    struct('name', 'Spiral',      'obj', ShapeSpiral(tp.Spiral),           'col', [0.9 0.55 0.1], 'T_final', 1*tp.T_duration_spiral,       'ground_size', 2.5), ...
    struct('name', 'Lissajous3D', 'obj', ShapeLissajous3D(tp.Lissajous3D), 'col', [0.2 0.7 0.3],  'T_final', 1*tp.T_duration_lissajous3D, 'ground_size', 3.0) ...
};

%% ODE OPTIONS
opts_ode = odeset('RelTol', 1e-7, 'AbsTol', 1e-8);

%% COMMON QUANTITIES
n_states  = 12;
T_hover = mp.m * mp.g;

%% MAIN LOOP OVER TRAJECTORIES
for ti = 1:numel(trajs)

    traj = trajs{ti};

    fprintf(' Trajectory: %s\n', traj.name);
    fprintf('Simulation time: %.3f s\n\n', traj.T_final);

    %% INITIAL CONDITION FROM FLATNESS MAP
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

    fprintf('Initial inputs (desired motor thrusts):\n');
    fprintf('  T = %.4f N   L = %.4f Nm   M = %.4f Nm   N = %.4f Nm\n\n', ...
        u0(1), u0(2), u0(3), u0(4));

    %% ODE INTEGRATION
    [t_sol, x_sol] = ode45( ...
        @(t, x) ode_dynamics(t, x, quad, traj.obj, mp, flatness_map), ...
        [0 traj.T_final], x0, opts_ode);

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

    %% FIGURE 1 — RIGID-BODY STATE TRACKING (states 1-12)
    state_names = {
        'x_I [m]','y_I [m]','z_I [m]',...
        'v_x [m/s]','v_y [m/s]','v_z [m/s]',...
        '\phi [deg]','\theta [deg]','\psi [deg]',...
        'p [rad/s]','q [rad/s]','r [rad/s]'};

    angle_scale = [1 1 1 1 1 1 180/pi 180/pi 180/pi 1 1 1];

    figure('Name', sprintf('%s State Tracking (A02)', traj.name), 'Color', 'w', 'Position', [30 50 1600 850]);
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
    sgtitle(sprintf('%s tracking (A02 Model)  |  max pos error = %.2e m', traj.name, max(pos_err)));


    %% FIGURE 2 — INPUTS (desired motor thrusts)
    figure('Name', sprintf('%s Inputs (A02)', traj.name), 'Color', 'w', 'Position', [30 50 900 500]);
    input_labels = {'T [N]', 'L [Nm]', 'M [Nm]', 'N [Nm]'};

    for j = 1:4
        subplot(4, 1, j);
        plot(t_sol, u_ref(:, j), 'Color', traj.col, 'LineWidth', 1.8); grid on;
        yline(T_hover, 'k--', 'T_{hover}');
        ylabel(input_labels{j});
    end
    xlabel('t [s]'); sgtitle(sprintf('%s reference inputs', traj.name));

    %% FIGURE 3 — POSITION ERROR
    figure('Name', sprintf('%s Position Error (A02)', traj.name), 'Color', 'w');
    plot(t_sol, pos_err, 'LineWidth', 2); grid on;
    xlabel('t [s]'); ylabel('||e_p|| [m]');
    title(sprintf('%s position tracking error (A02)', traj.name));

    %% FIGURE 4 — 3D DRONE ANIMATION
    fig3 = figure('Name', sprintf('3D %s Animation (A02)', traj.name), 'Color', 'w', 'Position', [900 100 800 750]);
    ax3 = axes(fig3);
    C00_PlotDrone.setup_axes(ax3, sprintf('3D Drone — %s Trajectory (A02)', traj.name));
    hold(ax3, 'on');

    sz = traj.ground_size;
    C00_PlotDrone.draw_ground(ax3, [-sz sz], [-sz sz], 0);

    plot3(ax3, x_ref(:, 1), x_ref(:, 2), x_ref(:, 3), 'k--', 'LineWidth', 1.5);
    h_traj = plot3(ax3, NaN, NaN, NaN, 'Color', traj.col, 'LineWidth', 1.5);
    h_drone = C00_PlotDrone.init(ax3, mp);

    lgd = findobj(fig3, 'Type', 'Legend');
    lgd.Location = 'southeast';

    pos0 = x_sol(1, 1:3)';
    R0 = Euler2R_ZYX(x_sol(1, 7), x_sol(1, 8), x_sol(1, 9));
    C00_PlotDrone.update(h_drone, pos0, R0);

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

        C00_PlotDrone.update(h_drone, pos_k, R_k);
        set(h_traj, 'XData', x_sol(1:k, 1), 'YData', x_sol(1:k, 2), 'ZData', x_sol(1:k, 3));

        title(ax3, sprintf([ ...
            '%s (C00)  |  t = %.2f s  |  ', ...
            '\\phi = %.1f°  ', ...
            '\\theta = %.1f°  ', ...
            '\\psi = %.1f°'], ...
            traj.name, t_sol(k), phi_k*180/pi, the_k*180/pi, psi_k*180/pi));
        drawnow;
    end

    fprintf('\n End simulation — %s \n\n', traj.name);

end

fprintf(' All trajectories done.\n');


%% LOCAL FUNCTIONS
function dxdt = ode_dynamics(t, x, quad, traj_obj, mp, flatness_map)
    % Reconstruct trajectory state and derivatives
    [s, ds, dds, ddds, dddds] = traj_obj.get_flat_outputs(t);

    % Get reference inputs (desired motor thrusts) from Flatness Map
    [~, uref] = flatness_map.map(mp, s, ds, dds, ddds, dddds);

    % Evaluate dynamics (ode)
    dxdt = quad.dynamics(t, x, uref);
end