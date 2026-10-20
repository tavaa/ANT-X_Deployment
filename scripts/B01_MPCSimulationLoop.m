% B01_MPCSimulationLoop.m
% LTV-MPC closed-loop simulation for the ANT-X quadrotor, B01 model.
% Scenarios: Circle, Lemniscate, Spiral
% Initial conditions: OnReference, Perturbation, Origin

clear; clc; close all;

%% PATHS
addpath(fullfile('models'));
addpath(fullfile('flatness'));
addpath(fullfile('configs'));
addpath(fullfile('control/MPC'));
addpath(fullfile('trajectories'));
addpath(fullfile('utils'));


%% OUTPUT DIRECTORY
[scriptDir, ~, ~] = fileparts(mfilename('fullpath'));
projectRoot    = fileparts(scriptDir);
resultsBaseDir = fullfile(projectRoot, 'results/MPC_B01/');
if ~exist(resultsBaseDir, 'dir'), mkdir(resultsBaseDir); end

todayStr = datestr(now, 'yyyy-mm-dd');
runID = 1;
while true
    runDirName = sprintf('%s_Run_%02d', todayStr, runID);
    fullRunDir = fullfile(resultsBaseDir, runDirName);
    if ~exist(fullRunDir, 'dir'), mkdir(fullRunDir); break; end
    runID = runID + 1;
end

fprintf('MPC Simulation Loop — ANT-X (B01, 16-state)\nOutput: %s\n\n', runDirName);


%% SETTINGS
pos_threshold = 0.05;   % [m] convergence detection threshold
hover_tail_duration = 3.0;   % [s]  Hover tail (Perturbation, Origin).
transmission_latency = 0.0005;   % [s] simulates delayed trasmission latency, not yet measured


%% SYSTEM
mp      = ModelParameters();
cp      = ControlParameters();
tp      = TrajectoryParameters();
Ts      = cp.Ts;
options = odeset('RelTol', 1e-8, 'AbsTol', 1e-9);

% model and controller
model = B01_ANTX_quadcopter(mp, true); % B01: single model instance, delay_motors=true
mpc   = B01_MPCController(model);

T_hover_pm = mp.m * mp.g / 4;   % [N] per-motor hover thrust

%% OSQP warm-up. 
% Throwaway generic hover-consistent problem, result discarded. 
x_warmup    = [zeros(12,1); T_hover_pm*ones(4,1)];
xref_warmup = zeros(12, cp.p + 2);
uref_warmup = T_hover_pm * ones(4, cp.p + 2);
mpc.init(x_warmup, xref_warmup, uref_warmup);
mpc.solve(x_warmup, xref_warmup(:, 1:cp.p+1), uref_warmup(:, 1:cp.p));


%% TRAJECTORIES (offline PWC reference generation, 12-state no-delay)
traj_circle      = ShapeCircle(tp.Circle);
traj_lemniscate2 = ShapeLemniscate2(tp.Lemniscate2);
traj_spiral      = ShapeSpiral(tp.Spiral);

[x_circle,      u_circle]      = B01_GeneratePWC_reference(traj_circle,      tp.T_duration_circle,      Ts, options);
[x_lemniscate2, u_lemniscate2] = B01_GeneratePWC_reference(traj_lemniscate2, tp.T_duration_lemniscate2, Ts, options);
[x_spiral,      u_spiral]      = B01_GeneratePWC_reference(traj_spiral,      tp.T_duration_spiral,      Ts, options);


%% SCENARIOS  { name, traj_obj, T_end, x_ref, u_ref }
scenarios = {
    'Circle',      traj_circle,      tp.T_duration_circle,      x_circle,      u_circle;
    'Lemniscate2', traj_lemniscate2, tp.T_duration_lemniscate2, x_lemniscate2, u_lemniscate2;
    'Spiral',      traj_spiral,      tp.T_duration_spiral,      x_spiral,      u_spiral;
};


%% JSON CONFIGURATION
cfg.model      = struct('name', class(model), 'flatness', 'B01_FlatnessMap', ...
                        'm', mp.m, 'g', mp.g, 'Jxx', mp.Jxx, 'Jyy', mp.Jyy, 'Jzz', mp.Jzz, ...
                        'Lp', mp.Lp, 'Mp', mp.Mq, 'Np', mp.Nr, ...
                        'tau_accel', mp.tau_accel, 'tau_decel', mp.tau_decel);
cfg.control    = struct('Ts', cp.Ts, 'horizon', cp.p, ...
                        'Q', cp.Q_diag_16, 'Qf', cp.Qf_diag_16, 'R', cp.R_diag_16, ...
                        'u_min', cp.u_min_16, 'u_max', cp.u_max_16, ...
                        'x_min', cp.x_min_16, 'x_max', cp.x_max_16);
cfg.trajectories = struct('Circle', tp.Circle, 'Lemniscate2', tp.Lemniscate2, 'Spiral', tp.Spiral);
fid = fopen(fullfile(fullRunDir, 'simulation_config.json'), 'w');
if fid ~= -1, fwrite(fid, jsonencode(cfg, 'PrettyPrint', true), 'char'); fclose(fid); end


%% LOG HEADERS
state_headers = {'step', ...
    'r_n','r_e','r_d','v_n','v_e','v_d','phi','theta','psi','p','q','r','T1','T2','T3','T4', ...
    'ref_rn','ref_re','ref_rd','ref_vn','ref_ve','ref_vd','ref_phi','ref_theta','ref_psi','ref_p','ref_q','ref_r','ref_T1','ref_T2','ref_T3','ref_T4'};

input_headers = {'step','T1','T2','T3','T4','T1_ref','T2_ref','T3_ref','T4_ref'};

error_headers = {'step', ...
    'position_error_m','attitude_error_deg','velocity_error_ms','rate_error_rads', ...
    'thrust_error_N', ...
    'yaw_error_raw_rad','yaw_error_corrected_rad','fullstate_error', ...
    'cumulative_mse_pos','cumulative_rmse_pos', ...
    'cumulative_mse_full','cumulative_rmse_full', ...
    'cumulative_mse_input','cumulative_rmse_input', ...
    'in_track','solve_time_s'};

state_labels = {'$r_n$ [m]','$r_e$ [m]','$r_d$ [m]', ...
    '$v_n$ [m/s]','$v_e$ [m/s]','$v_d$ [m/s]', ...
    '$\phi$ [deg]','$\theta$ [deg]','$\psi$ [deg]', ...
    '$p$ [rad/s]','$q$ [rad/s]','$r$ [rad/s]', ...
    '$T_1$ [N]','$T_2$ [N]','$T_3$ [N]','$T_4$ [N]'};

input_labels = {'T_1^d [N]','T_2^d [N]','T_3^d [N]','T_4^d [N]'};

angle_idx = [7 8 9];


%% MAIN LOOP
initial_conditions = {'OnReference','Perturbation','Origin'};
%initial_conditions = {'Origin'};
metrics_summary  = {};
metrics_extended = {};

for c = 1:numel(initial_conditions)
    cond_name = initial_conditions{c};

    for s = 1:size(scenarios, 1)

        scenario_name = scenarios{s, 1};
        traj_obj      = scenarios{s, 2};
        T_end_base    = scenarios{s, 3};   % end of the REFERENCE 
        x_ref_pwc     = scenarios{s, 4};   % [12 x T_ref], no-delay reference
        u_ref_pwc     = scenarios{s, 5};   % [4  x T_ref]

        run_dir = fullfile(fullRunDir, scenario_name, cond_name);
        if ~exist(run_dir, 'dir'), mkdir(run_dir); end

        fprintf('%s | %s\n', scenario_name, cond_name);

        %% Hover at the end of the reference perturbation / origin scenario.
        if strcmp(cond_name, 'Perturbation') || strcmp(cond_name, 'Origin')
            n_tail = round(hover_tail_duration / Ts);
            x_hover_pt = x_ref_pwc(:, end);
            x_hover_pt(4:6)   = 0;
            x_hover_pt(7:8)   = 0;
            x_hover_pt(10:12) = 0;
            x_ref_use = [x_ref_pwc, repmat(x_hover_pt, 1, n_tail)];
            u_ref_use = [u_ref_pwc, repmat(T_hover_pm*ones(4,1), 1, n_tail)];
            T_end_use = T_end_base + hover_tail_duration;
        else
            x_ref_use = x_ref_pwc;
            u_ref_use = u_ref_pwc;
            T_end_use = T_end_base;
        end

        %% Initial condition
        [s0,ds0,dds0,ddds0,dddds0] = traj_obj.get_flat_outputs(0);
        [x0_nom_16, ~] = B01_FlatnessMap.map(mp, s0, ds0, dds0, ddds0, dddds0);
        current_x = x0_nom_16(:);

        if strcmp(cond_name, 'OnReference')
            mpc.init(current_x, x_ref_use, u_ref_use);

        elseif strcmp(cond_name, 'Perturbation')
            current_x(1)     = current_x(1) + 0.20;
            current_x(2)     = current_x(2) - 0.20;
            current_x(4:6)   = 0;
            current_x(7:8)   = 0;
            current_x(9)     = 0;
            current_x(10:12) = 0;
            mpc.init(current_x, x_ref_use, u_ref_use);

        elseif strcmp(cond_name, 'Origin')
            current_x        = zeros(16, 1);
            current_x(3)     = -0.02;
            current_x(13:16) = cp.u_min_16;
            if strcmp(scenario_name, 'Lemniscate2')
                current_x(1) = 0.75;
            end
            mpc.init(current_x, x_ref_use, u_ref_use);
        end

        x_hover      = x_ref_use(:, end);
        x_hover(4:6) = 0; x_hover(7:8) = 0; x_hover(10:12) = 0;
        u_hover      = T_hover_pm * ones(4,1);

        t_steps      = 0 : Ts : T_end_use;
        N_steps      = length(t_steps);
        N_steps_base = length(0 : Ts : T_end_base);  

        log_states = zeros(N_steps, numel(state_headers));
        log_inputs = zeros(N_steps, numel(input_headers));
        log_errors = zeros(N_steps, numel(error_headers));
        log_perst  = zeros(N_steps, 16);

        sum_sq_pos  = 0;
        sum_sq_full = 0;
        sum_sq_u    = 0;
        solve_times = zeros(N_steps, 1);

        u_prev = T_hover_pm * ones(4,1);


        %% RECEDING HORIZON
        for k = 1:N_steps
            t_now = t_steps(k);

            x_meas = current_x;

            xref_seq = zeros(12, cp.p + 1);
            uref_seq = zeros(4,  cp.p);
            for j = 0:cp.p
                idx = k + j;
                if idx <= size(x_ref_use, 2)
                    xref_seq(:, j+1) = x_ref_use(:, idx);
                    if j < cp.p, uref_seq(:, j+1) = u_ref_use(:, idx); end
                else
                    xref_seq(:, j+1) = x_hover;
                    if j < cp.p, uref_seq(:, j+1) = u_hover; end
                end
            end

            xref_seq_qp = wrap_yaw_reference(xref_seq, x_meas(9));

            tic;
            [u_opt, ~] = mpc.solve(x_meas, xref_seq_qp, uref_seq);
            calc_time  = toc;
            solve_times(k) = calc_time;
            
            % total actuation_delay
            total_actuation_delay = calc_time + transmission_latency;

            if total_actuation_delay < Ts
                [~, x_ode1] = ode45(@(t,x) model.dynamics(t, x, u_prev), ...
                    [t_now, t_now + total_actuation_delay], current_x, options);
                x_mid = x_ode1(end, :)';
                [~, x_ode2] = ode45(@(t,x) model.dynamics(t, x, u_opt), ...
                    [t_now + total_actuation_delay, t_now + Ts], x_mid, options);
                next_x = x_ode2(end, :)';
            else
                warning('MPC:ActuationOverrun', ...
                    'total_actuation_delay=%.1fms exceeds Ts=%.1fms at t=%.2fs -- u_opt skipped this step.', ...
                    total_actuation_delay*1e3, Ts*1e3, t_now);
                [~, x_ode1] = ode45(@(t,x) model.dynamics(t, x, u_prev), ...
                    [t_now, t_now + Ts], current_x, options);
                next_x = x_ode1(end, :)';
            end

            u_prev = u_opt;

            ref_now_12 = xref_seq(:, 1);
            uref_now   = uref_seq(:, 1);
            ref_now_16 = [ref_now_12; uref_now];

            e = current_x - ref_now_16;
            e_yaw_raw  = e(9);
            e_yaw_wrap = atan2(sin(e_yaw_raw), cos(e_yaw_raw));
            e(9) = e_yaw_wrap;

            err_pos    = norm(e(1:3));
            err_vel    = norm(e(4:6));
            err_att    = norm(e(7:9)) * 180/pi;
            err_rate   = norm(e(10:12));
            err_thrust = norm(e(13:16));
            err_full   = norm(e);
            err_u      = norm(u_opt - uref_now);

            sum_sq_pos  = sum_sq_pos  + err_pos^2;
            sum_sq_full = sum_sq_full + err_full^2;
            sum_sq_u    = sum_sq_u    + err_u^2;

            in_track = double(err_pos < pos_threshold);

            log_states(k,:) = [k, current_x', ref_now_16'];
            log_inputs(k,:) = [k, u_opt', uref_now'];
            log_errors(k,:) = [k, err_pos, err_att, err_vel, err_rate, err_thrust, ...
                e_yaw_raw, e_yaw_wrap, err_full, ...
                sum_sq_pos/k, sqrt(sum_sq_pos/k), ...
                sum_sq_full/k, sqrt(sum_sq_full/k), ...
                sum_sq_u/k, sqrt(sum_sq_u/k), ...
                in_track, calc_time];

            e_ps = abs(e);
            e_ps(angle_idx) = e_ps(angle_idx) * 180/pi;
            log_perst(k,:) = e_ps';

            current_x = next_x;
        end


        %% POST-LOOP ANALYSIS
        in_track_vec = logical(log_errors(1:N_steps_base, 16));

        k_conv = NaN;
        for ki = 1:N_steps_base-4
            if all(in_track_vec(ki:ki+4)), k_conv = ki; break; end
        end
        t_conv = NaN;
        if ~isnan(k_conv), t_conv = t_steps(k_conv); end

        if sum(in_track_vec) > 10
            rmse_ss_pos  = sqrt(mean(log_errors(in_track_vec, 2).^2));
            rmse_ss_full = sqrt(mean(log_errors(in_track_vec, 9).^2));
        else
            rmse_ss_pos = NaN; rmse_ss_full = NaN;
        end

        max_pos_err = max(log_errors(1:N_steps_base, 2));
        max_yaw_err = max(abs(log_errors(1:N_steps_base, 7)));
        max_att_err = max(log_errors(1:N_steps_base, 3));
        avg_solve   = mean(solve_times) * 1e3;   
        worst_solve = max(solve_times)  * 1e3;
        std_solve   = std(solve_times)  * 1e3;
        n_slow      = sum(solve_times > Ts * 0.8);

        perstate_rmse = zeros(1,16);
        for si = 1:16, perstate_rmse(si) = sqrt(mean(log_perst(1:N_steps_base, si).^2)); end

        % Reference-only final RMSE/MSE
        rmse_pos_ref   = log_errors(N_steps_base, 11);
        rmse_full_ref  = log_errors(N_steps_base, 13);
        rmse_input_ref = log_errors(N_steps_base, 15);
        mse_pos_ref    = log_errors(N_steps_base, 10);
        mse_full_ref   = log_errors(N_steps_base, 12);
        mse_input_ref  = log_errors(N_steps_base, 14);

        fprintf('  steps=%d (ref=%d)  avg=%.1fms  worst=%.1fms  slow=%d\n', ...
            N_steps, N_steps_base, avg_solve, worst_solve, n_slow);
        fprintf('  RMSE_pos=%.4fm  RMSE_ss=%.4fm  t_conv=%.2fs \n\n', ...
            rmse_pos_ref, rmse_ss_pos, t_conv);


        %% FIGURES
        fig_title = sprintf('%s — %s', scenario_name, cond_name);

        %% Figure 1: State Tracking -- FULL duration (tail included)
        fig1 = figure('Name', sprintf('States: %s', fig_title), ...
            'Color','w', 'Position',[50 50 1400 950]);
        for i = 1:16
            subplot(4,4,i);
            sc = 1; if ismember(i, angle_idx), sc = 180/pi; end
            plot(t_steps, log_states(:,i+1)*sc,  'b',   'LineWidth',1.5); hold on; grid on;
            plot(t_steps, log_states(:,i+17)*sc, 'r--', 'LineWidth',1.2);
            xline(T_end_base, 'k:', 'LineWidth', 0.8);
            ylabel(state_labels{i}, 'Interpreter','latex');
            xlim([0 T_end_use]);
            if i >= 13, xlabel('t [s]'); end
            if i == 1, legend('Simulated','Reference','Location','best','Interpreter','none'); end
        end
        sgtitle(sprintf('States Tracking: %s', fig_title), 'Interpreter','none');
        saveas(fig1, fullfile(run_dir, 'plot_states.png'));


        %% Figure 2: Applied Inputs -- FULL duration (tail included)
        fig2 = figure('Name', sprintf('Inputs: %s', fig_title), ...
            'Color','w', 'Position',[150 100 800 600]);
        u_ylims = [cp.u_min_16(1), cp.u_max_16(1);
                   cp.u_min_16(2), cp.u_max_16(2);
                   cp.u_min_16(3), cp.u_max_16(3);
                   cp.u_min_16(4), cp.u_max_16(4)];
        for i = 1:4
            subplot(2,2,i);
            plot(t_steps, log_inputs(:,i+1), 'b',   'LineWidth',1.5); hold on; grid on;
            plot(t_steps, log_inputs(:,i+5), 'r--', 'LineWidth',1.2);
            xline(T_end_base, 'k:', 'LineWidth', 0.8);
            yline(u_ylims(i,1), 'k:', 'LineWidth',0.8);
            yline(u_ylims(i,2), 'k:', 'LineWidth',0.8);
            yline(T_hover_pm,   'g:', 'LineWidth',0.8);
            ylim(u_ylims(i,:));
            ylabel(input_labels{i}); xlim([0 T_end_use]);
            if i >= 3, xlabel('t [s]'); end
            if i == 1, legend('Applied','Reference','Location','best','Interpreter','none'); end
        end
        sgtitle(sprintf('Inputs Applied: %s', fig_title), 'Interpreter','none');
        saveas(fig2, fullfile(run_dir, 'plot_inputs.png'));


        %% Figure 3: Tracking Analysis -- REFERENCE SPAN ONLY (no tail)
        t_steps_ref  = t_steps(1:N_steps_base);
        solve_ref    = solve_times(1:N_steps_base);

        fig3 = figure('Name', sprintf('Analysis: %s', fig_title), ...
            'Color','w', 'Position',[250 100 1000 900]);

        subplot(3,2,1);
        semilogy(t_steps_ref, log_errors(1:N_steps_base,2)+1e-12, 'b', 'LineWidth',1.5); hold on; grid on;
        yline(pos_threshold, 'r--', sprintf('%.0f cm', pos_threshold*100), ...
            'LabelVerticalAlignment','bottom');
        if ~isnan(t_conv)
            xline(t_conv, 'k-', sprintf('t_{conv} = %.1fs', t_conv), ...
                'LabelOrientation','horizontal');
        end
        ylabel('Position error [m]'); xlabel('t [s]');
        title('Position Error (reference span)');

        subplot(3,2,2);
        plot(t_steps_ref, log_errors(1:N_steps_base,3), 'Color',[0.8 0.2 0.2], 'LineWidth',1.5); grid on;
        ylabel('Attitude error [deg]'); xlabel('t [s]');
        title('Attitude Error  (\phi, \theta, \psi)');

        subplot(3,2,3);
        plot(t_steps_ref, log_errors(1:N_steps_base,7)*180/pi, 'Color',[0.6 0.6 0.6], 'LineWidth',1.0); hold on; grid on;
        plot(t_steps_ref, log_errors(1:N_steps_base,8)*180/pi, 'Color',[0.1 0.4 0.8], 'LineWidth',1.8);
        legend('algebraic','corrected','Location','best');
        ylabel('[deg]'); xlabel('t [s]');
        title('Yaw Error');

        subplot(3,2,4);
        plot(t_steps_ref, log_errors(1:N_steps_base,6), 'Color',[0.5 0.2 0.6], 'LineWidth',1.5); grid on;
        ylabel('Thrust error [N]'); xlabel('t [s]');
        title('Motor Thrust Error  (T_1..T_4, delay-induced)');

        subplot(3,2,5);
        c_pos    = [0.25 0.45 0.80];
        c_vel    = [0.20 0.65 0.30];
        c_ang    = [0.85 0.50 0.10];
        c_rate   = [0.75 0.20 0.20];
        c_thrust = [0.55 0.25 0.65];
        hold on; grid on;
        h1 = bar(1:3,   perstate_rmse(1:3),   'FaceColor', c_pos);
        h2 = bar(4:6,   perstate_rmse(4:6),   'FaceColor', c_vel);
        h3 = bar(7:9,   perstate_rmse(7:9),   'FaceColor', c_ang);
        h4 = bar(10:12, perstate_rmse(10:12), 'FaceColor', c_rate);
        h5 = bar(13:16, perstate_rmse(13:16), 'FaceColor', c_thrust);
        xticks(1:16);
        ax5 = gca;
        ax5.TickLabelInterpreter = 'latex';
        ax5.XTickLabel = { ...
            '$r_n$', '$r_e$', '$r_d$', ...
            '$v_n$', '$v_e$', '$v_d$', ...
            '$\phi$', '$\theta$', '$\psi$', ...
            '$p$', '$q$', '$r$', ...
            '$T_1$','$T_2$','$T_3$','$T_4$'};
        xtickangle(30);
        ylabel('RMSE'); xlabel('State');
        title('Per-State RMSE (reference span)');
        legend([h1 h2 h3 h4 h5], {'position','velocity','angles','rates','thrust'}, ...
            'Location','northeast','FontSize',7);

        subplot(3,2,6);
        plot(t_steps_ref, solve_ref*1e3, 'Color',[0.3 0.3 0.3], 'LineWidth',0.8); hold on; grid on;
        yline(Ts*1e3,     'r--', sprintf('T_s = %.0f ms', Ts*1e3), ...
            'LabelHorizontalAlignment','left');
        yline(Ts*0.8*1e3, 'r:',  '80% budget', ...
            'LabelHorizontalAlignment','right');
        ylabel('Solve time [ms]'); xlabel('t [s]');
        title(sprintf('QP Solve Time  (worst = %.1f ms, ref span)', max(solve_ref)*1e3));

        sgtitle(sprintf('Tracking Analysis: %s (reference span only)', fig_title), 'Interpreter','none');
        saveas(fig3, fullfile(run_dir, 'plot_analysis.png'));


        %% GIF
        gif_name  = fullfile(run_dir, 'tracking.gif');
        GenerateTrackingSimulationGIF(log_states, x_ref_use, gif_name, fig_title, '2D', mp, PlotDrone2);


        %% CSV -- full duration logs kept complete 
        writetable(array2table(log_states, 'VariableNames', state_headers), ...
            fullfile(run_dir, 'states.csv'));
        writetable(array2table(log_inputs, 'VariableNames', input_headers), ...
            fullfile(run_dir, 'inputs.csv'));
        writetable(array2table(log_errors, 'VariableNames', error_headers), ...
            fullfile(run_dir, 'errors.csv'));
        ps_names = arrayfun(@(i) sprintf('state_%d_error', i), 1:16, 'UniformOutput', false);
        writetable(array2table(log_perst, 'VariableNames', ps_names), ...
            fullfile(run_dir, 'per_state_errors.csv'));


        %% ACCUMULATE METRICS -- reference span only
        metrics_summary(end+1,:) = {cond_name, scenario_name, ...
            rmse_pos_ref, rmse_full_ref, rmse_input_ref, ...
            avg_solve, N_steps_base, T_end_base};

        metrics_extended(end+1,:) = {cond_name, scenario_name, ...
            rmse_pos_ref, rmse_ss_pos, ...
            rmse_full_ref, rmse_ss_full, ...
            mse_pos_ref, mse_full_ref, ...
            t_conv, ...
            max_pos_err, max_yaw_err*180/pi, max_att_err, ...
            rmse_input_ref, mse_input_ref, ...
            avg_solve, worst_solve, std_solve, n_slow, ...
            N_steps_base, T_end_base};

    end
end


%% EXPORT SUMMARY
cols_sum = {'initial_condition','scenario', ...
    'RMSE_position_m','RMSE_fullstate','RMSE_input', ...
    'avg_solve_time_ms','steps_reference_span','reference_duration_s'};
writetable(cell2table(metrics_summary, 'VariableNames', cols_sum), ...
    fullfile(fullRunDir, 'MPC_metrics_summary.csv'));


%% EXPORT EXTENDED
cols_ext = {'initial_condition','scenario', ...
    'rmse_position_all_steps_m','rmse_position_steadystate_m', ...
    'rmse_fullstate_all_steps','rmse_fullstate_steadystate', ...
    'mse_position_m2','mse_fullstate', ...
    'convergence_time_s', ...
    'peak_position_error_m','peak_yaw_error_deg','peak_attitude_error_deg', ...
    'rmse_input','mse_input', ...
    'avg_solve_time_ms','worst_solve_time_ms','std_solve_time_ms', ...
    'steps_exceeding_80pct_budget','steps_reference_span','reference_duration_s'};
writetable(cell2table(metrics_extended, 'VariableNames', cols_ext), ...
    fullfile(fullRunDir, 'MPC_metrics_extended.csv'));

fprintf('Done. Results: %s\n', fullRunDir);

%% DEBUG: Omega_i plausibility check across all 9 (IC x scenario) combinations
% Converts logged T1..T4 [N] to Omega1..Omega4 [rad/s] via Omega=sqrt(T/K_T).

K_T = mp.K_T;
Omega_min = sqrt(cp.u_min_16(1) / K_T);
Omega_max = sqrt(cp.u_max_16(1) / K_T);

motor_colors = {[0.20 0.45 0.80], [0.85 0.20 0.20], [0.20 0.65 0.30], [0.60 0.30 0.70]};

figure('Name', 'Omega_i plausibility check (all combos)', 'Color', 'w', 'Position', [50 50 1500 1100]);
tl = tiledlayout(3, 3, 'Padding', 'compact', 'TileSpacing', 'compact');

for c = 1:numel(initial_conditions)
    cond_name = initial_conditions{c};
    for s = 1:size(scenarios, 1)
        scenario_name = scenarios{s, 1};

        run_dir  = fullfile(fullRunDir, scenario_name, cond_name);
        Tstates  = readtable(fullfile(run_dir, 'states.csv'));
        t_       = (Tstates.step - 1) * Ts;

        Om1 = sqrt(max(Tstates.T1, 0) / K_T);
        Om2 = sqrt(max(Tstates.T2, 0) / K_T);
        Om3 = sqrt(max(Tstates.T3, 0) / K_T);
        Om4 = sqrt(max(Tstates.T4, 0) / K_T);

        nexttile;
        plot(t_, Om1, 'Color', motor_colors{1}, 'LineWidth', 1.1); hold on; grid on;
        plot(t_, Om2, 'Color', motor_colors{2}, 'LineWidth', 1.1);
        plot(t_, Om3, 'Color', motor_colors{3}, 'LineWidth', 1.1);
        plot(t_, Om4, 'Color', motor_colors{4}, 'LineWidth', 1.1);
        yline(Omega_min, 'k:', 'LineWidth', 0.7);
        yline(Omega_max, 'k:', 'LineWidth', 0.7);

        title(sprintf('%s — %s', scenario_name, cond_name), 'FontSize', 9);
        if s == 1, ylabel('\Omega [rad/s]'); end
        if c == numel(initial_conditions), xlabel('t [s]'); end
    end
end

lgd = legend({'\Omega_1','\Omega_2','\Omega_3','\Omega_4'}, 'Orientation', 'horizontal');
lgd.Layout.Tile = 'south';
sgtitle(tl, 'Motor angular velocities \Omega_i = sqrt(T_i / K_T)');