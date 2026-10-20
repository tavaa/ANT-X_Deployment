%% A02_EXPORT_REFERENCE
% Generate the NLP reference (A02 model) for the four shapes and export it
% to results/reference_A02/<shape>.csv
%
% Columns: t, xI, yI, zI, vx, vy, vz, phi, theta, psi, p, q, r, T, L, M, N.
% Rows: nodes t_k = k*Ts, k = 0..N (last input repeated).

clear; clc; close all;

%% PATHS
addpath(fullfile('models'));
addpath(fullfile('flatness'));
addpath(fullfile('configs'));
addpath(fullfile('trajectories'));
addpath(fullfile('utils'));
addpath(fullfile('casadi'));
addpath(fullfile('deployment'));

%% SETUP
mp = ModelParameters();
cp = ControlParameters();
tp = TrajectoryParameters();
Ts = cp.Ts;

out_dir = fullfile('results', 'reference_A02');
if ~exist(out_dir, 'dir'), mkdir(out_dir); end

% name, trajectory object, duration, periodic
shapes = {
    'Circle',     ShapeCircle(tp.Circle),            tp.T_duration_circle,       true
    'Lemniscate', ShapeLemniscate2(tp.Lemniscate),   tp.T_duration_lemniscate,   true
    'Spiral',     ShapeSpiral(tp.Spiral),            tp.T_duration_spiral,       false
    'Lissajous',  ShapeLissajous3D(tp.Lissajous3D),  tp.T_duration_lissajous3D,  true
};

var_names = {'t', 'xI','yI','zI', 'vx','vy','vz', 'phi','theta','psi', 'p','q','r', 'T','L','M','N'};

%% GENERATE, VALIDATE, EXPORT
for i = 1:size(shapes, 1)
    name     = shapes{i, 1};
    traj_obj = shapes{i, 2};
    T_end    = shapes{i, 3};
    periodic = shapes{i, 4};

    % create optimized trajectory
    nlp_opts = struct('periodic', periodic, 'print_level', 0);
    [x_ref, u_ref, t_steps, info] = A02_GeneratePWC_reference(traj_obj, T_end, Ts, nlp_opts);
    
    % validate trajectory
    rep = A02_Validate_Trajectory(x_ref, u_ref, t_steps, info, name, false);
    
    % not export if not validated correctly
    if ~rep.pass
        warning('A02_Export_Reference:ValidationFail', '%s: validation FAILED, not exported.', name);
        continue;
    end

    data = [t_steps(:), x_ref.', u_ref.'];
    writetable(array2table(data, 'VariableNames', var_names), fullfile(out_dir, [name '.csv']));

    fprintf('%-10s | N = %d | Ts = %.2f s | %s -> %s.csv\n', ...
        name, numel(t_steps) - 1, Ts, info.stats.return_status, name);
end