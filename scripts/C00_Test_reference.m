%% Test optimized PWC reference. C00 model, all trajectories, Ts = 0.1 s.
%
%  Produce, validate and export the offline reference that the MPC will
%  track. For every trajectory listed below the script
%    1. solves the NLP of C00_GeneratePWC_reference (u_bar + du*),
%    2. validates the result with C00_TrajectoryValidation
%       (thrust limits, ode45 replay on the fine grid, MPC state box,
%        tracking of x^r(t) within the Bryson tolerances),
%    3. exports it as CSV,
%    4. plots it against the CONTINUOUS reference of C00_FlatnessMap.
%  A summary table at the end shows which trajectory passes and which
%  check fails.
%
%  EXPORT (ground truth = the N+1 nodes of the NLP, NO RK4 sub-intervals)
%      results/reference/trajectory_<name>.csv   one file per trajectory
%  Columns: t, r_n r_e r_d v_n v_e v_d phi theta psi p q r T1 T2 T3 T4, u1..u4
%  (SI units; u_i = desired motor thrusts; the input of the LAST row repeats
%  u_{N-1} as padding and is never applied).

clear; clc; close all;

%% PATHS
addpath(fullfile('models'));
addpath(fullfile('flatness'));
addpath(fullfile('configs'));
addpath(fullfile('trajectories'));
addpath(fullfile('utils'));
addpath(fullfile('casadi'));

%% SETTINGS
mp = ModelParameters;
tp = TrajectoryParameters;

Ts         = 0.10;                            % GCS deployment sampling time [s]
export_only_if_valid = false;                
out_dir    = fullfile('results', 'reference');
T_hover_pm = mp.m * mp.g / 4;

%% TRAJECTORY LIST
% periodic = closed curve (x_N = x^ref(t_N) is imposed).
trajs = { ...
    makeEntry('Circle',      ShapeCircle(tp.Circle),           tp.T_duration_circle,      true,  [0.10 0.45 0.90]), ...
    makeEntry('Lemniscate2', ShapeLemniscate2(tp.Lemniscate2), tp.T_duration_lemniscate2, true,  [0.60 0.20 0.80]), ...
    makeEntry('Spiral',      ShapeSpiral(tp.Spiral),           tp.T_duration_spiral,      false, [0.90 0.55 0.10]), ...
    makeEntry('Lissajous3D', ShapeLissajous3D(tp.Lissajous3D), tp.T_duration_lissajous3D, true,  [0.10 0.70 0.35]) ...
    };

fprintf('C00 optimized reference test | Ts = %.3f s\n\n', Ts);

%% MAIN LOOP
summary = cell(numel(trajs), 1);

for ti = 1:numel(trajs)

    traj = trajs{ti};

    fprintf(' Trajectory: %s | T = %.2f s (%d steps) | periodic = %d\n', ...
        traj.name, traj.T_final, round(traj.T_final / Ts), traj.periodic);

    %% SOLVE THE NLP
    tic;
    [x_ref, u_ref, t_steps, info] = C00_GeneratePWC_reference( ...
        traj.obj, traj.T_final, Ts, struct('periodic', traj.periodic));
    t_solve = toc;

    %% VALIDATE
    report = C00_TrajectoryValidation(x_ref, u_ref, t_steps, info, traj.name);

    %% EXPORT 
    saved = '-';
    if report.pass || ~export_only_if_valid
        file_traj = fullfile(out_dir, sprintf('trajectory_%s.csv', lower(traj.name)));
        saveReferenceCSV(file_traj, t_steps, x_ref, u_ref);
        saved = file_traj;
        fprintf(' Saved: %s\n', file_traj);
        if ~report.pass
            warning('C00_Test_reference:SavedDespiteFail', ...
                '%s: validation FAILED but the CSV was exported anyway.', traj.name);
        end
    else
        warning('C00_Test_reference:NotSaved', ...
            '%s: validation FAILED -> NOT exported (export_only_if_valid = true).', traj.name);
    end

    %% PLOTS
    plotReference(traj, t_steps, x_ref, u_ref, info, T_hover_pm);

    %% Summary row
    failed = report.checks(strcmp(report.checks(:, 4), 'FAIL'), 1);
    summary{ti} = struct( ...
        'name',    traj.name, ...
        'verdict', ternary(report.pass, 'PASS', 'FAIL'), ...
        'n_warn',  report.n_warn, ...
        'solve_s', t_solve, ...
        'max_du',  1e3 * max(abs(info.du(:))), ...
        'headroom',report.headroom.dn_min, ...
        'saved',   saved, ...
        'failed',  {failed});
end

%% SUMMARY TABLE

fprintf(' SUMMARY (Ts = %.3f s)\n', Ts);
fprintf(' %-13s %-7s %-5s %-9s %-11s %-13s %s\n', ...
    'trajectory', 'verdict', 'warn', 'solve [s]', 'max|du| mN', 'min room [N]', 'saved');
for ti = 1:numel(summary)
    s = summary{ti};
    fprintf(' %-13s %-7s %-5d %-9.1f %-11.2f %-13.3f %s\n', ...
        s.name, s.verdict, s.n_warn, s.solve_s, s.max_du, s.headroom, s.saved);
    for f = 1:numel(s.failed)
        fprintf('     -> FAILED: %s\n', s.failed{f});
    end
end


%% LOCAL FUNCTIONS
function e = makeEntry(name, obj, T_final, periodic, col)
% MAKEENTRY  One row of the trajectory list.
    e = struct('name', name, 'obj', obj, 'T_final', T_final, 'periodic', periodic, 'col', col);
end

function out = ternary(cond, a, b)
% TERNARY  Inline if.
    if cond, out = a; else, out = b; end
end

function plotReference(traj, t_steps, x_ref, u_ref, info, T_hover_pm)
% PLOTREFERENCE  Optimized nodes vs CONTINUOUS flatness reference:
%                position, 3D path, inputs (u_bar vs u_bar + du*).

    t_dense  = info.t_dense;
    xr_dense = info.xr_dense;

    % Position
    figure('Name', sprintf('%s Position (C00 NLP)', traj.name), 'Color', 'w', 'Position', [30 50 900 650]);
    pos_labels = {'r_n [m]', 'r_e [m]', 'r_d [m]'};
    for j = 1:3
        subplot(3, 1, j);
        plot(t_dense, xr_dense(j, :), 'k-', 'LineWidth', 1.2); hold on;
        plot(t_steps, x_ref(j, :), 'o-', 'Color', traj.col, 'LineWidth', 1.4, 'MarkerSize', 3);
        ylabel(pos_labels{j}); grid on;
        if j == 1, legend('continuous ref (flatness)', 'optimized nodes', 'Location', 'best'); end
    end
    xlabel('t [s]');
    sgtitle(sprintf('%s: optimized reference vs continuous flatness reference', traj.name));

    % 3D path
    figure('Name', sprintf('%s 3D Path (C00 NLP)', traj.name), 'Color', 'w', 'Position', [950 50 700 650]);
    plot3(xr_dense(1, :), xr_dense(2, :), xr_dense(3, :), 'k-', 'LineWidth', 1.2); hold on;
    plot3(x_ref(1, :), x_ref(2, :), x_ref(3, :), 'o-', 'Color', traj.col, 'LineWidth', 1.4, 'MarkerSize', 3);
    set(gca, 'ZDir', 'reverse', 'YDir', 'reverse');
    grid on; axis equal;
    xlabel('r_n [m]'); ylabel('r_e [m]'); zlabel('r_d [m]');
    legend('continuous ref', 'optimized nodes', 'Location', 'best');
    title(sprintf('%s: 3D path (NED)', traj.name));

    % Inputs
    figure('Name', sprintf('%s Inputs (C00 NLP)', traj.name), 'Color', 'w', 'Position', [30 50 900 600]);
    input_labels = {'T_1 [N]', 'T_2 [N]', 'T_3 [N]', 'T_4 [N]'};
    for j = 1:4
        subplot(4, 1, j);
        stairs(t_steps, info.u_bar(j, :), 'k--', 'LineWidth', 1.2); hold on;
        stairs(t_steps, u_ref(j, :), 'Color', traj.col, 'LineWidth', 1.6);
        yline(T_hover_pm, 'k:', 'T_{hover}/4');
        ylabel(input_labels{j}); grid on;
        if j == 1, legend('u_{bar} (case B)', 'u_{bar}+\Delta u^* (case C)', 'Location', 'best'); end
    end
    xlabel('t [s]');
    sgtitle(sprintf('%s: interval-averaged input vs NLP-corrected input', traj.name));
end

function saveReferenceCSV(file, t_steps, x_ref, u_ref)
% SAVEREFERENCECSV  Writes the NLP reference (nodes only) as CSV.
%
%   One row per node k = 0..N: time, state (nx columns), input (4 columns).
%   The input of the last row repeats u_{N-1} (padding, never applied).
%   Full double precision (%.12e), header row with column names.

    folder = fileparts(file);
    if ~isempty(folder) && ~exist(folder, 'dir')
        mkdir(folder);
    end

    nx = size(x_ref, 1);
    state_names = {'r_n','r_e','r_d','v_n','v_e','v_d','phi','theta','psi', ...
                   'p','q','r','T1','T2','T3','T4'};
    names = [{'t'}, state_names(1:nx), {'u1','u2','u3','u4'}];

    D = [t_steps(:), x_ref.', u_ref.'];          % (N+1) x (1 + nx + 4)

    fid = fopen(file, 'w');
    if fid < 0
        error('saveReferenceCSV:Open', 'Cannot open %s for writing.', file);
    end
    fprintf(fid, '%s\n', strjoin(names, ','));
    fprintf(fid, [repmat('%.12e,', 1, size(D, 2) - 1), '%.12e\n'], D.');
    fclose(fid);
end