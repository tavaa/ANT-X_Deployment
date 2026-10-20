%% OPENLOOP_PWC_REFERENCE
%
% Compares the open-loop trajectory obtained with:
%   - Continuous feedforward inputs  u(t)  -> ode45 with time-varying u
%   - Piece-Wise Constant inputs     u_k   -> u held constant at each k*Ts
%
% For each scenario (Circle, Lemniscate, Spiral) and for each Ts value
% the script:
%   - Runs the continuous reference simulation
%   - Calls GeneratePWC_reference for every Ts
%   - Plots trajectory comparison (2D or 3D) in a 2x2 subplot figure

clear; clc; close all;

%% PATHS
addpath(fullfile('models'));
addpath(fullfile('flatness'));
addpath(fullfile('configs'));
addpath(fullfile('trajectories'));
addpath(fullfile('utils'));

%% SETTINGS
use_quaternions = false; % Set to true for A01 (quaternion), false for A00 (Euler Z-Y-X)

%% PARAMETERS
%Ts_values = [0.001, 0.01, 0.02, 0.05];   % sampling periods to test
Ts_values = [0.01, 0.02, 0.05, 0.1];   % sampling periods to test
colors    = {'b',  'g',  'm',  'r'};  

% Model parameters struct
mp = ModelParameters();

% Dynamic model and flatness map selection
if use_quaternions
    model    = A01_SimplifiedModel_quat(mp);
    flatness = A01_FlatnessMap();
    fprintf('Selected model: A01 (13-State Quaternion representation)\n');
else
    model    = A00_SimplifiedModel(mp);
    flatness = A00_FlatnessMap();
    fprintf('Selected model: A00 (12-State Euler Z-Y-X representation)\n');
end

% ode45 tolerances 
options = odeset('RelTol', 1e-8, 'AbsTol', 1e-9);

%% TRAJECTORY CONFIGURATIONS
tp = TrajectoryParameters();

traj_circle      = ShapeCircle(tp.Circle);
traj_lemniscate  = ShapeLemniscate(tp.Lemniscate);
traj_lemniscate2 = ShapeLemniscate2(tp.Lemniscate2);
traj_spiral      = ShapeSpiral(tp.Spiral);

%% SCENARIOS 
scenarios = {
    'Circle',      traj_circle,      tp.T_duration_circle,      '2D';
    'Lemniscate',  traj_lemniscate,  tp.T_duration_lemniscate,  '2D';
    'Lemniscate2', traj_lemniscate2, tp.T_duration_lemniscate2, '2D';
    'Spiral',      traj_spiral,      tp.T_duration_spiral,      '3D';
};

fprintf("Starting Continuous vs PWC Inputs model simulation.\n")

%% MAIN LOOP OVER SCENARIOS
for s = 1 : size(scenarios, 1)

    shape_name = scenarios{s, 1};
    traj_obj   = scenarios{s, 2};
    T_sim      = scenarios{s, 3};
    plot_type  = scenarios{s, 4};

    fprintf('\nSimulating Scenario: %s\n', shape_name);

    % CONTINUOUS REFERENCE SIMULATION
    fprintf('Continuous simulation\n');

    % Initial condition from flatness map at t = 0
    [s0, ds0, dds0, ddds0, dddds0] = traj_obj.get_flat_outputs(0);
    [x0, ~] = flatness.map(mp, s0, ds0, dds0, ddds0, dddds0);

    % Integrate with time-varying continuous input
    dyn_cont = @(t, x) continuous_dynamics(t, x, model, flatness, traj_obj, mp, use_quaternions);
    [t_cont, x_cont] = ode45(dyn_cont, [0, T_sim], x0, options);

    % PWC REFERENCE SIMULATIONS 
    pwc_results = cell(1, length(Ts_values));

    for i = 1 : length(Ts_values)
        curr_Ts = Ts_values(i);
        fprintf('PWC Inputs simulation  Ts = %.3f s \n', curr_Ts);

        [x_ref, u_ref, t_steps] = GeneratePWC_reference( ...
            traj_obj, T_sim, curr_Ts, model, flatness, options);

        pwc_results{i}.Ts    = curr_Ts;
        pwc_results{i}.t     = t_steps;
        pwc_results{i}.x     = x_ref;
        pwc_results{i}.u     = u_ref;
    end
    
    %% PLOTTING — 2x2 subplots, one per Ts
    figure('Name',     sprintf('%s — Continuous vs PWC', shape_name), ...
           'Color',    'w', ...
           'Position', [100, 100, 1200, 900]);

    for i = 1 : length(Ts_values)
        res = pwc_results{i};
        ax  = subplot(2, 2, i);
        set(ax, 'ZDir', 'reverse');
        set(ax, 'YDir', 'reverse');
        hold(ax, 'on');
        grid(ax, 'on');
        box(ax,  'on');

        if strcmp(plot_type, '2D')
            % continuous reference (black dashed)
            plot(ax, x_cont(:, 1), x_cont(:, 2), ...
                 'k--', 'LineWidth', 1.2, 'DisplayName', 'Continuous');
            % PWC trajectory 
            plot(ax, res.x(1, :), res.x(2, :), ...
                 '-', 'Color', colors{i}, 'LineWidth', 2.0, ...
                 'DisplayName', sprintf('PWC Inputs'));
            xlabel(ax, 'X [m]');
            ylabel(ax, 'Y [m]');

        else   % '3D'
            plot3(ax, x_cont(:, 1), x_cont(:, 2), x_cont(:, 3), ...
                  'k--', 'LineWidth', 1.2, 'DisplayName', 'Continuous');
            plot3(ax, res.x(1, :), res.x(2, :), res.x(3, :), ...
                  '-', 'Color', colors{i}, 'LineWidth', 2.0, ...
                  'DisplayName', sprintf('PWC Inputs'));
            xlabel(ax, 'X [m]');
            ylabel(ax, 'Y [m]');
            zlabel(ax, 'Z [m]');
            view(ax, 30, 30);
        end

        axis(ax, 'equal');
        title(ax,  sprintf('Ts = %.3f s', res.Ts));
        legend(ax, 'Location', 'best');
    end

    sgtitle(sprintf('%s — Open-Loop Divergence: Continuous vs PWC Inputs', shape_name), ...
            'FontSize', 13, 'FontWeight', 'bold');

end

fprintf('\nDone.\n');

%% LOCAL FUNCTION — continuous dynamics 
function dxdt = continuous_dynamics(t, x, model, flatness, traj_obj, mp, use_quaternions)
    % Flat outputs and their derivatives at time t
    [s0, ds0, dds0, ddds0, dddds0] = traj_obj.get_flat_outputs(t);

    % Feedforward input from flatness map
    [~, u_cont] = flatness.map(mp, s0, ds0, dds0, ddds0, dddds0);

    % Normalize quaternion to prevent numerical integration drift
    if use_quaternions
        q = x(7:10);
        x(7:10) = q / norm(q);
    end

    % Nonlinear model dynamics
    dxdt = model.dynamics(t, x, u_cont);
end