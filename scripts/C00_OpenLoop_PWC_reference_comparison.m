%% OPEN-LOOP COMPARISON OF PWC REFERENCE-GENERATION STRATEGIES
%
% Compare three PWC input strategies on the same ANT-X model:
%
%   A : u_k = u^r(t_k)
%   B : u_k = interval-average input
%   C : u_k = interval-average input + NLP correction
%
% All inputs are applied open loop from x^r(0) and compared with the
% continuous reference trajectory.
%
% Strategies:
%   A  continuous flatness input sampled at the start of each interval
%   B  ODE-based interval average
%   C  NLP-corrected PWC input
%
% The plant is the 12-state padded model with motor delay disabled.
%
% For each shape, sampling time and number of laps:
%   1. Generate the NLP PWC reference.
%   2. Build inputs A, B and C.
%   3. Propagate all three inputs with ode45.
%   4. Compute position tracking errors.
%   5. Run trajectory validation.
%
% Results are saved in results/comparison/ and reused through a cache.
%
% Output figures:
%   snapshot_<shape>       final trajectory comparison
%   tracking_error_<shape> position error over time
%   summary                maximum position error versus Ts
%
% Ts_list and laps_list are kept unchanged.

clear; clc; close all;

%% PATHS
addpath(fullfile('models'));
addpath(fullfile('flatness'));
addpath(fullfile('configs'));
addpath(fullfile('trajectories'));
addpath(fullfile('utils'));
addpath(fullfile('casadi'));

%% SETTINGS
shapes_list = {'Lemniscate', 'Lissajous3D'};   % Lemniscate = Bernoulli (tp.Lemniscate)
Ts_list     = [0.01, 0.02, 0.05, 0.10];         % [s]  
laps_list   = [1, 2, 3];                        % 

sig_p   = 0.01;          % [m] position tolerance
tau_c   = 0.20;          % [s] NLP correction time scale
h_nlp   = 0.01;          % [s] target RK4 step inside the NLP

div_pos  = 3.0;          % [m] stop if position error exceeds this value
div_tilt = deg2rad(60);  % [rad] stop if roll or pitch exceeds this value

snap_tile_w     = 300;   % [px] snapshot tile width
snap_tile_h_max = 190;   % [px] maximum snapshot tile height

use_cache = true;        % reuse cached results
check_B01 = true;        % compare strategy B with B01 once

opts_ode = odeset('RelTol', 1e-9, 'AbsTol', 1e-10);
out_root = fullfile('results', 'comparison');

%% STRATEGY STYLE
S.tag  = {'A', 'B', 'C'};
S.name = {'A: $u^r(t_k)$ held over $T_s$', ...
          'B: $\bar u_k$ (interval mean)', ...
          'C: $\bar u_k + \Delta u_k^\star$ (NLP)'};
S.col  = [0.85 0.33 0.10;     % A vermillion
          0.93 0.69 0.13;     % B amber
          0.00 0.45 0.74];    % C blue

%% SETUP
mp    = ModelParameters;
tp    = TrajectoryParameters;
model = C00_ANTX_quadcopter(mp, false);   % 12-state plant, no motor lag
Ts_max = max(Ts_list);
rows   = {};
dirs   = makeDirs(out_root);

fprintf('PWC reference comparison (padded, no-lag plant) | shapes: %s | Ts: %s s | laps: %s\n\n', ...
    strjoin(shapes_list, ', '), mat2str(Ts_list), mat2str(laps_list));

E_max_all = nan(3, numel(Ts_list), numel(laps_list), numel(shapes_list));

%% MAIN LOOP -- SHAPE -> Ts -> LAPS
for si = 1:numel(shapes_list)
    [traj_obj, T_lap, omega] = makeShape(shapes_list{si}, tp, Ts_max);

    case_grid = cell(numel(Ts_list), numel(laps_list));

    for ti = 1:numel(Ts_list)
        Ts = Ts_list(ti);
        M  = max(2, round(Ts / h_nlp));

        for li = 1:numel(laps_list)
            L     = laps_list(li);
            T_sim = L * T_lap;
            key   = sprintf('%s_Ts%03d_L%d', shapes_list{si}, round(1000*Ts), L);
            cfg   = [sig_p, tau_c, M, omega, T_sim, div_pos, div_tilt];
            cache_file = fullfile(dirs.cache, [key '.mat']);

            %% Run or load the case
            R = [];
            if use_cache && isfile(cache_file)
                Cc = load(cache_file, 'R');
                if isfield(Cc.R, 'cfg') && isequal(Cc.R.cfg, cfg)
                    R = Cc.R;
                end
            end
            if isempty(R)
                R = runCase(traj_obj, T_sim, T_lap, L, Ts, M, model, mp, ...
                            sig_p, tau_c, opts_ode, div_pos, div_tilt, key, S);
                R.cfg = cfg;
                save(cache_file, 'R');
            end

            fprintf(['%-11s Ts=%.2f L=%d | max|e_p| [m]  A %8.2e  B %8.2e  C %8.2e', ...
                     ' | validation  A %-8s B %-8s C %-8s | NLP %6.1f s\n'], ...
                shapes_list{si}, Ts, L, R.e_max(1), R.e_max(2), R.e_max(3), ...
                R.verdict{1}, R.verdict{2}, R.verdict{3}, R.t_nlp);

            E_max_all(:, ti, li, si) = R.e_max;
            case_grid{ti, li} = R;
            rows(end+1, :) = {shapes_list{si}, Ts, L, R.e_max, R.verdict, R.t_nlp}; %#ok<SAGROW>

            %% B01 consistency check
            if check_B01 && Ts == Ts_max && L == laps_list(1)
                checkB01(traj_obj, T_sim, Ts, M, opts_ode, R);
            end
        end
    end

    %% Grid figures for this shape
    plotSnapshotGrid(case_grid, Ts_list, laps_list, S, shapes_list{si}, ...
                      fullfile(dirs.png, ['snapshot_' lower(shapes_list{si})]), ...
                      snap_tile_w, snap_tile_h_max);
    plotErrorGrid(case_grid, Ts_list, laps_list, S, shapes_list{si}, sig_p, ...
                  fullfile(dirs.png, ['tracking_error_' lower(shapes_list{si})]));
end

%% SUMMARY
plotSummary(E_max_all, Ts_list, laps_list, shapes_list, S, ...
            fullfile(dirs.png, 'summary'), sig_p);
printTable(rows);


%% LOCAL FUNCTIONS

function R = runCase(traj_obj, T_sim, T_lap, L, Ts, M, model, mp, sig_p, tau_c, ...
                     opts_ode, div_pos, div_tilt, key, S)
% RUNCASE  Generate, simulate and validate strategies A, B and C.
%
% A: flatness input sampled at each interval start.
% B: interval-average input.
% C: NLP-corrected PWC input.
%
% All three strategies use the same plant and ode45 propagation.

    %% Generate NLP reference and interval-average input
    nlp_opts = struct('delay_motors', false, 'M', M, 'periodic', true, ...
                      'sig_p', sig_p, 'tau_c', tau_c, 'print_level', 0);
    x_C = []; u_C = []; t_steps = []; info = []; %#ok<NASGU>
    tic;
    evalc('[x_C, u_C, t_steps, info] = C00_GeneratePWC_reference(traj_obj, T_sim, Ts, nlp_opts);');
    t_nlp = toc;

    N  = numel(t_steps) - 1;
    x0 = info.xr_dense(:, 1);

    %% Build input sequences
    U = cell(1, 3);

    % A: continuous reference input at the start of each interval
    U{1} = zeros(4, N);
    for k = 1:N
        [s, ds, dds, ddds, dddds] = traj_obj.get_flat_outputs(t_steps(k));
        [~, uk] = C00_FlatnessMap.map(mp, s, ds, dds, ddds, dddds);
        U{1}(:, k) = uk;
    end

    % B: interval-average input
    U{2} = info.u_bar(:, 1:N);

    % C: NLP-corrected input
    U{3} = u_C(:, 1:N);

    %% Open-loop propagation
    P = cell(1, 3);  Xn = cell(1, 3);  t_div = nan(1, 3);
    for s = 1:3
        [Xd, Xn{s}, t_div(s)] = propagatePWC(model, U{s}, x0, t_steps, M, opts_ode, ...
                                             info.xr_dense, div_pos, div_tilt);
        P{s} = Xd(1:3, :);
    end

    pr = info.xr_dense(1:3, :);

    % Position error with respect to the continuous reference
    E  = nan(3, size(pr, 2));
    for s = 1:3
        E(s, :) = vecnorm(P{s} - pr, 2, 1);
    end
    e_max = max(E, [], 2, 'omitnan');

    % Error at the end of each lap
    e_lap = nan(3, L);
    for l = 1:L
        idx = round(l * T_lap / Ts) * M + 1;
        e_lap(:, l) = E(:, min(idx, size(E, 2)));
    end

    %% Validation
    verdict = cell(1, 3);  n_fail = nan(1, 3);  ratio_pos = nan(1, 3);

    for s = 1:3
        if ~isnan(t_div(s))
            verdict{s} = 'DIVERGED';
            continue;
        end

        info_s = info;

        if s < 3
            info_s.du         = U{s} - U{2};
            info_s.stats      = struct('success', true, 'return_status', 'n/a (forward sim)');
            info_s.max_defect = 0;
            x_s = Xn{s};
        else
            x_s = x_C;
        end

        u_s  = [U{s}, U{s}(:, end)];
        name = sprintf('%s-%s', key, S.tag{s}); %#ok<NASGU>
        rep = [];

        try
            evalc('rep = C00_TrajectoryValidation(x_s, u_s, t_steps, info_s, name, false);');

            if rep.pass
                verdict{s} = 'PASS';
            else
                verdict{s} = 'FAIL';
            end

            n_fail(s)    = rep.n_fail;
            ratio_pos(s) = rep.ratio.position;
        catch
            verdict{s} = 'ERROR';
        end
    end

    %% Pack results
    R.t_steps    = t_steps;
    R.t_dense    = info.t_dense;
    R.pr         = pr;
    R.P          = P;
    R.E          = E;
    R.t_div      = t_div;
    R.e_max      = e_max;
    R.e_lap      = e_lap;
    R.verdict    = verdict;
    R.n_fail     = n_fail;
    R.ratio_pos  = ratio_pos;
    R.t_nlp      = t_nlp;
    R.nlp_status = info.stats.return_status;
    R.du_max     = max(abs(info.du(:)));
    R.T_sim      = T_sim;
    R.T_lap      = T_lap;
    R.laps       = L;
    R.Ts         = Ts;
    R.M          = M;
end


function [Xd, Xn, t_div] = propagatePWC(model, U, x0, t_steps, M, opts_ode, xr_dense, div_pos, div_tilt)
% PROPAGATEPWC  Propagate a PWC input with ode45.
%
% U(:,k) is constant over [t_k, t_k+Ts].
% Xd contains the fine-grid trajectory.
% Xn contains the states at the PWC nodes.
% The simulation stops if the state becomes invalid or diverges.

    N  = size(U, 2);
    nx = numel(x0);
    h  = (t_steps(2) - t_steps(1)) / M;

    Xd = nan(nx, N*M + 1);
    Xd(:, 1) = x0;
    x = x0;
    t_div = NaN;

    w_state = warning('off', 'MATLAB:ode45:IntegrationTolNotMet');

    for k = 1:N
        tk   = t_steps(k) + (0:M) * h;
        cols = (k-1)*M + 1 : k*M + 1;

        try
            [~, xs] = ode45(@(t, x) model.dynamics(t, x, U(:, k)), tk, x, opts_ode);
        catch
            t_div = tk(1);
            break;
        end

        if size(xs, 1) ~= M + 1
            t_div = tk(1);
            break;
        end

        Xd(:, cols) = xs.';
        x = xs(end, :).';

        e_pos = norm(x(1:3) - xr_dense(1:3, cols(end)));

        if ~all(isfinite(x)) || e_pos > div_pos || max(abs(x(7:8))) > div_tilt
            t_div = tk(end);
            break;
        end
    end

    warning(w_state);

    Xn = Xd(:, 1:M:end);
end


function [obj, T_lap, omega] = makeShape(name, tp, Ts_max)
% MAKESHAPE  Create the requested trajectory and lap time.

    switch name
        case 'Lemniscate'
            p = tp.Lemniscate;
            obj = ShapeLemniscate(p);

        case 'Lissajous3D'
            p = tp.Lissajous3D;
            obj = ShapeLissajous3D(p);

        otherwise
            error('C00_OpenLoop_PWC_reference_comparison:Shape', 'Unknown shape "%s".', name);
    end

    omega = p.omega;

    T_lap = Ts_max * round((2*pi / omega) / Ts_max);

    if abs(T_lap - 2*pi/omega) > 1e-6
        fprintf(['  note: %s lap 2*pi/omega = %.4f s used as T_lap = %.2f s ', ...
                 '(set omega = 2*pi/%.2f for exact closure)\n'], name, 2*pi/omega, T_lap);
    end
end


function dirs = makeDirs(out_root)
% MAKEDIRS  Create output folders.

    dirs = struct('png', fullfile(out_root, 'png'), 'cache', fullfile(out_root, 'cache'));

    f = fieldnames(dirs);
    for i = 1:numel(f)
        if ~exist(dirs.(f{i}), 'dir')
            mkdir(dirs.(f{i}));
        end
    end
end


function checkB01(traj_obj, T_sim, Ts, M, opts_ode, R)
% CHECKB01  Check that strategy B matches the B01 generator.

    try
        [x_B01, ~, ~] = B01_GeneratePWC_reference(traj_obj, T_sim, Ts, opts_ode);

        pB = R.P{2}(:, 1:M:end);
        n  = min(size(x_B01, 2), size(pB, 2));

        d  = max(vecnorm(x_B01(1:3, 1:n) - pB(:, 1:n), 2, 1), [], 'omitnan');

        fprintf('  [B01 check] B01_GeneratePWC_reference vs strategy B: max |dp| = %.2e m\n', d);

    catch err
        fprintf('  [B01 check] skipped: %s\n', err.message);
    end
end


function plotSnapshotGrid(case_grid, Ts_list, laps_list, S, shape_name, file_base, tile_w, tile_h_max)
% PLOTSNAPSHOTGRID  Plot all open-loop trajectories in one grid.
%
% Rows = Ts, columns = number of laps.
% The continuous reference is drawn first.
% A, B and C are dashed with different styles.
% 3D shapes use a 3D NED view; planar shapes use a top view.

    nT = numel(Ts_list);
    nL = numel(laps_list);

    G  = snapGeometry(case_grid{1, 1}.pr, tile_w, tile_h_max);

    sty = {'--', '-.', '-'};
    lw  = [2.2, 2.0, 1.8];
    
    fig = figure('Color', 'w', 'Position', ...
                 [60 40, round(nL*G.tile_w + 110), round(nT*(G.tile_h + 16) + 135)]);

    tl  = tiledlayout(fig, nT, nL, 'Padding', 'compact', 'TileSpacing', 'compact');

    axs = gobjects(nT, nL);
    h_leg = gobjects(1, 4);

    for ti = 1:nT
        for li = 1:nL
            R  = case_grid{ti, li};

            ax = nexttile(tl, (ti-1)*nL + li);
            axs(ti, li) = ax;

            hold(ax, 'on');
            grid(ax, 'on');
            box(ax, 'on');
            ax.GridAlpha = 0.18;
            ax.MinorGridAlpha = 0.08;
            ax.LineWidth = 0.8;

            % Continuous reference
            if G.use3d
                h_ref = plot3(ax, R.pr(1, :), R.pr(2, :), R.pr(3, :), ...
                              'k-', 'LineWidth', 0.8);
            else
                h_ref = plot3(ax, R.pr(1, :), R.pr(2, :), R.pr(3, :), ...
                'Color', [0.10 0.10 0.10], 'LineStyle', '-', 'LineWidth', 1.8);
            end

            % Strategies A, B and C
            h_tr = gobjects(1, 3);

            for s = 1:3
                P = R.P{s};

                if G.use3d
                    h_tr(s) = plot3(ax, P(1, :), P(2, :), P(3, :), sty{s}, ...
                                    'Color', S.col(s, :), 'LineWidth', lw(s));
                else
                    h_tr(s) = plot(ax, P(G.horiz, :), P(G.vert, :), sty{s}, ...
                                   'Color', S.col(s, :), 'LineWidth', lw(s));
                end

                % Mark the point where a run diverged
                j = find(~isnan(P(1, :)), 1, 'last');

                if ~isnan(R.t_div(s)) && ~isempty(j)
                    if G.use3d
                        plot3(ax, P(1, j), P(2, j), P(3, j), 'x', ...
                              'Color', S.col(s, :), 'MarkerSize', 8, 'LineWidth', 2);
                    else
                        plot(ax, P(G.horiz, j), P(G.vert, j), 'x', ...
                             'Color', S.col(s, :), 'MarkerSize', 8, 'LineWidth', 2);
                    end
                end

                if G.use3d
                    plot3(ax, P(1,1), P(2,1), P(3,1), 'o', ...
                     'Color', S.col(s,:), ...
                     'MarkerFaceColor', S.col(s,:), ...
                     'MarkerSize', 4, 'HandleVisibility', 'off');

                    plot3(ax, P(1,end), P(2,end), P(3,end), 's', ...
                     'Color', S.col(s,:), ...
                     'MarkerFaceColor', 'w', ...
                     'MarkerSize', 4, 'HandleVisibility', 'off');
                end
            end

            if ti == 1 && li == 1
                h_leg = [h_ref, h_tr];
            end

            % Axis limits from the reference
            if G.use3d
                xlim(ax, [G.lo(1), G.hi(1)]);
                ylim(ax, [G.lo(2), G.hi(2)]);
                zlim(ax, [G.lo(3), G.hi(3)]);

                set(ax, 'ZDir', 'reverse', 'YDir', 'reverse');
                view(ax, [-42 28]);
                camproj(ax, 'perspective');
                camlight(ax, 'headlight');
                lighting(ax, 'gouraud');

            else
                xlim(ax, [G.lo(G.horiz), G.hi(G.horiz)]);
                ylim(ax, [G.lo(G.vert), G.hi(G.vert)]);

                if G.horiz == 1
                    set(ax, 'YDir', 'reverse');
                end
            end

            daspect(ax, [1 1 1]);

            % Show labels only on the outer tiles
            if G.use3d
                if ti == nT && li == 1
                    xlabel(ax, 'r_n [m]');
                    ylabel(ax, 'r_e [m]');
                    zlabel(ax, 'r_d [m]');
                else
                    set(ax, 'XTickLabel', [], 'YTickLabel', [], 'ZTickLabel', []);
                end
            else
                if ti < nT
                    set(ax, 'XTickLabel', []);
                else
                    xlabel(ax, G.lab{G.horiz});
                end

                if li > 1
                    set(ax, 'YTickLabel', []);
                else
                    ylabel(ax, G.lab{G.vert});
                end
            end

            title(ax, sprintf('T_s = %.2f s, %d lap(s)', Ts_list(ti), laps_list(li)), ...
                  'FontSize', 7, 'FontWeight', 'normal');

            set(ax, 'FontSize', 6);
        end
    end

    lgd = legend(axs(1, 1), h_leg, [{'continuous reference $x^r(t)$'}, S.name], ...
                 'Interpreter', 'latex', 'Orientation', 'horizontal', 'NumColumns', 2, ...
                 'FontSize', 8);

    lgd.ItemTokenSize = [30 18];
    lgd.Layout.Tile = 'south';

    title(tl, {sprintf('%s ,  final-state snapshot', shape_name)}, 'FontSize', 9);

    exportgraphics(fig, [file_base '.png'], 'Resolution', 200);
    exportgraphics(fig, [file_base '.pdf'], 'ContentType', 'vector');
end


function G = snapGeometry(pr, tile_w, tile_h_max)
% SNAPGEOMETRY  Choose 2D/3D view and plot limits from the reference.

    pad = 0.15;
    ext = max(pr, [], 2) - min(pr, [], 2);

    G.lab   = {'r_n [m]', 'r_e [m]', 'r_d [m]'};
    G.use3d = ext(3) > 0.05;
    G.lo    = min(pr, [], 2) - pad;
    G.hi    = max(pr, [], 2) + pad;

    if G.use3d
        G.horiz = 1;
        G.vert = 2;
        G.a = 1.3;

    else
        if ext(2) > ext(1)
            G.horiz = 2;
            G.vert = 1;
        else
            G.horiz = 1;
            G.vert = 2;
        end

        G.a = (ext(G.horiz) + 2*pad) / (ext(G.vert) + 2*pad);
    end

    G.tile_w = tile_w;
    G.tile_h = tile_w / G.a;

    if G.tile_h > tile_h_max
        G.tile_h = tile_h_max;
        G.tile_w = G.a * tile_h_max;
    end
end


function plotErrorGrid(case_grid, Ts_list, laps_list, S, shape_name, sig_p, file_base)
% PLOTERRORGRID  Plot position error for all cases.
%
% Rows = Ts, columns = number of laps.
% The y-axis is logarithmic and sigma_p is shown as a reference line.

    nT = numel(Ts_list);
    nL = numel(laps_list);

    fig = figure('Color', 'w', 'Position', [80 60 260*nL + 120, 220*nT + 140]);

    tl  = tiledlayout(fig, nT, nL, 'Padding', 'compact', 'TileSpacing', 'compact');

    axs = gobjects(nT, nL);

    for ti = 1:nT
        for li = 1:nL
            R  = case_grid{ti, li};

            ax = nexttile(tl, (ti-1)*nL + li);
            axs(ti, li) = ax;

            hold(ax, 'on');
            grid(ax, 'on');
            box(ax, 'on');

            set(ax, 'YScale', 'log');

            h_e = gobjects(1, 3);

            for s = 1:3
                h_e(s) = plot(ax, R.t_dense, R.E(s, :), '-', ...
                              'Color', S.col(s, :), 'LineWidth', 1.3);
            end

            yline(ax, sig_p, 'k--');

            for l = 1:R.laps - 1
                xline(ax, l * R.T_lap, 'k:');
            end

            xlim(ax, [0, R.T_sim]);
            ylim(ax, [1e-9, 1e1]);

            if li == 1
                ylabel(ax, sprintf('T_s = %.2f s\n||e_p|| [m]', Ts_list(ti)));
            else
                ylabel(ax, '');
            end

            if ti == nT
                xlabel(ax, 't [s]');
            end

            if ti == 1
                title(ax, sprintf('%d lap(s)', laps_list(li)));
            end

            set(ax, 'FontSize', 7);
        end
    end

    lgd = legend(axs(1, 1), h_e, S.name, ...
                 'Interpreter', 'latex', 'Orientation', 'horizontal', 'FontSize', 8);

    lgd.Layout.Tile = 'south';

    title(tl, sprintf('%s ,  tracking error vs continuous reference.', shape_name));

    exportgraphics(fig, [file_base '.png'], 'Resolution', 200);
    exportgraphics(fig, [file_base '.pdf'], 'ContentType', 'vector');
end


function plotSummary(E_max, Ts_list, laps_list, shapes_list, S, file_base, sig_p)
% PLOTSUMMARY  Plot maximum position error versus Ts.
%
% Colors identify strategies.
% Line styles identify the number of laps.

    styles = {'-o', '--s', ':^', '-.d'};

    fig = figure('Color', 'w', 'Position', [100 100 1150 450]);

    tl  = tiledlayout(fig, 1, numel(shapes_list), ...
                      'Padding', 'compact', 'TileSpacing', 'compact');

    for si = 1:numel(shapes_list)
        ax = nexttile(tl);
        hold(ax, 'on');
        grid(ax, 'on');
        box(ax, 'on');

        set(ax, 'XScale', 'log', 'YScale', 'log');

        for s = 1:3
            for li = 1:numel(laps_list)
                plot(ax, Ts_list, squeeze(E_max(s, :, li, si)), ...
                     styles{min(li, numel(styles))}, ...
                     'Color', S.col(s, :), ...
                     'MarkerFaceColor', S.col(s, :), ...
                     'LineWidth', 1.4);
            end
        end

        yline(ax, sig_p, 'k--', '\sigma_p');

        xticks(ax, Ts_list);
        xlim(ax, [0.8*min(Ts_list), 1.25*max(Ts_list)]);

        xlabel(ax, 'T_s [s]');
        ylabel(ax, 'max ||e_p|| [m]');
        title(ax, shapes_list{si});
    end

    % Legend: colors = strategies, line styles = laps
    h = gobjects(1, 3 + numel(laps_list));

    for s = 1:3
        h(s) = plot(ax, NaN, NaN, '-', ...
                    'Color', S.col(s, :), 'LineWidth', 2);
    end

    names = S.name;

    for li = 1:numel(laps_list)
        h(3 + li) = plot(ax, NaN, NaN, ...
                         styles{min(li, numel(styles))}, ...
                         'Color', 'k', 'LineWidth', 1.2);

        names{end+1} = sprintf('%d lap(s)', laps_list(li)); %#ok<AGROW>
    end

    lgd = legend(ax, h, names, ...
                 'Interpreter', 'latex', ...
                 'Orientation', 'horizontal', ...
                 'NumColumns', 3 + numel(laps_list));

    lgd.Layout.Tile = 'south';

    title(tl, 'Open-loop max position error vs T_s');

    exportgraphics(fig, [file_base '.png'], 'Resolution', 200);
    exportgraphics(fig, [file_base '.pdf'], 'ContentType', 'vector');
end


function printTable(rows)
% PRINTTABLE  Print a compact summary of all cases.


    fprintf(' %-11s %-5s %-4s | %-9s %-9s %-9s | %-8s %-8s %-8s | %s\n', ...
        'shape', 'Ts', 'laps', 'e_A [m]', 'e_B [m]', 'e_C [m]', ...
        'A', 'B', 'C', 'NLP [s]');

    fprintf('------------------------------------------------------------------------------------------\n');

    for i = 1:size(rows, 1)
        r = rows(i, :);

        fprintf(' %-11s %-5.2f %-4d | %-9.2e %-9.2e %-9.2e | %-8s %-8s %-8s | %.1f\n', ...
            r{1}, r{2}, r{3}, ...
            r{4}(1), r{4}(2), r{4}(3), ...
            r{5}{1}, r{5}{2}, r{5}{3}, r{6});
    end

end