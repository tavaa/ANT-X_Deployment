% MPCSimulationLoop_13states.m
% LTV-MPC simulation for the ANT-X quadrotor.
% Scenarios: Circle, Lemniscate, Spiral
% Initial conditions: OnReference, Perturbation, Origin
% Quaternions simulation

clear; clc; close all;

%% PATHS
addpath(fullfile('models'));
addpath(fullfile('flatness'));
addpath(fullfile('configs'));
addpath(fullfile('control/MPC'));
addpath(fullfile('trajectories'));
addpath(fullfile('utils'));


%% OUTPUT DIRECTORY
[scriptDir, ~, ...
    ~] = fileparts(mfilename('fullpath'));
projectRoot    = fileparts(scriptDir);
resultsBaseDir = fullfile(projectRoot, 'results/MPC/');
if ~exist(resultsBaseDir, 'dir'), mkdir(resultsBaseDir); end

todayStr = datestr(now, 'yyyy-mm-dd');
runID = 1;
while true
    runDirName = sprintf('%s_Run_%02d', todayStr, runID);
    fullRunDir = fullfile(resultsBaseDir, runDirName);
    if ~exist(fullRunDir, 'dir'), mkdir(fullRunDir); break; end
    runID = runID + 1;
end

fprintf('MPC Simulation Loop 13 states — ANT-X\nOutput: %s\n\n', runDirName);


%% SETTINGS
pos_threshold = 0.05;   % [m] convergence detection threshold


%% SYSTEM
mp      = ModelParameters();
cp      = ControlParameters();
tp      = TrajectoryParameters();
Ts      = cp.Ts;
options = odeset('RelTol', 1e-8, 'AbsTol', 1e-9);
model    = A01_SimplifiedModel_quat(mp);
flatness = A01_FlatnessMap();
mpc      = MPCController_13states(model);


%% TRAJECTORIES (offline PWC reference generation)
traj_circle      = ShapeCircle(tp.Circle);
traj_lemniscate2 = ShapeLemniscate2(tp.Lemniscate2);
traj_spiral      = ShapeSpiral(tp.Spiral);

% T_extra ensures the MPC horizon never exceeds the reference length.
T_extra = cp.p * Ts + Ts;

[x_circle,      u_circle]      = GeneratePWC_reference(traj_circle,      tp.T_duration_circle      + T_extra, Ts, model, flatness, options);
[x_lemniscate2, u_lemniscate2] = GeneratePWC_reference(traj_lemniscate2, tp.T_duration_lemniscate2 + T_extra, Ts, model, flatness, options);
[x_spiral,      u_spiral]      = GeneratePWC_reference(traj_spiral,      tp.T_duration_spiral      + T_extra, Ts, model, flatness, options);


%% SCENARIOS  { name, traj_obj, T_end, x_ref, u_ref }
scenarios = {
    'Circle',      traj_circle,      tp.T_duration_circle,      x_circle,      u_circle;
    'Lemniscate2', traj_lemniscate2, tp.T_duration_lemniscate2, x_lemniscate2, u_lemniscate2;
    'Spiral',      traj_spiral,      tp.T_duration_spiral,      x_spiral,      u_spiral;
};


%% JSON CONFIGURATION
cfg.model      = struct('name', class(model), 'flatness', class(flatness), ...
                        'm', mp.m, 'g', mp.g, 'Jxx', mp.Jxx, 'Jyy', mp.Jyy, 'Jzz', mp.Jzz);
cfg.control    = struct('Ts', cp.Ts, 'horizon', cp.p, ...
                        'Q', cp.Q_diag_13, 'Qf', cp.Qf_diag_13, 'R', cp.R_diag, ...
                        'u_min', cp.u_min, 'u_max', cp.u_max, ...
                        'x_min', cp.x_min_13, 'x_max', cp.x_max_13);
cfg.trajectories = struct('Circle', tp.Circle, 'Lemniscate2', tp.Lemniscate2, 'Spiral', tp.Spiral);
fid = fopen(fullfile(fullRunDir, 'simulation_config.json'), 'w');
if fid ~= -1, fwrite(fid, jsonencode(cfg, 'PrettyPrint', true), 'char'); fclose(fid); end


%% LOG HEADERS
state_headers = {'step', ...
    'x_I','y_I','z_I','u_B','v_B','w_B','q_0','q_1','q_2', 'q_3', 'p','q','r', ...
    'ref_x','ref_y','ref_z','ref_xb','ref_yb','ref_zb','ref_q0','ref_q1','ref_q2', 'ref_q3', 'ref_p','ref_q','ref_r'};

input_headers = {'step','T','L','M','N','T_ref','L_ref','M_ref','N_ref'};

error_headers = {'step', ...
    'position_error_m','attitude_error_deg','velocity_error_ms','rate_error_rads', ...
    'yaw_error_raw_rad','yaw_error_corrected_rad','fullstate_error', ...
    'cumulative_mse_pos','cumulative_rmse_pos', ...
    'cumulative_mse_full','cumulative_rmse_full', ...
    'cumulative_mse_input','cumulative_rmse_input', ...
    'in_track','solve_time_s'};

state_labels = {'$x_I$ [m]','$y_I$ [m]','$z_I$ [m]', ...
    '$\dot{x}_B$ [m/s]','$\dot{y}_B$ [m/s]','$\dot{z}_B$ [m/s]', ...
    '$q_0$','$q_1$','$q_2$', '$q_3$', ...
    '$p$ [rad/s]','$q_r$ [rad/s]','$r$ [rad/s]'};

input_labels = {'T — Thrust [N]','L — Roll [Nm]','M — Pitch [Nm]','N — Yaw [Nm]'};


%% MAIN LOOP
initial_conditions = {'OnReference','Perturbation','Origin'};
metrics_summary  = {};
metrics_extended = {};

for c = 1:numel(initial_conditions)
    cond_name = initial_conditions{c};

    for s = 1:size(scenarios, 1)

        scenario_name = scenarios{s, 1};
        traj_obj      = scenarios{s, 2};
        T_end         = scenarios{s, 3};
        x_ref_pwc     = scenarios{s, 4};
        u_ref_pwc     = scenarios{s, 5};

        run_dir = fullfile(fullRunDir, scenario_name, cond_name);
        if ~exist(run_dir, 'dir'), mkdir(run_dir); end

        fprintf('%s | %s\n', scenario_name, cond_name);

        % Initial condition
        [s0,ds0,dds0,ddds0,dddds0] = traj_obj.get_flat_outputs(0);
        [x0_nom, ~] = flatness.map(mp, s0, ds0, dds0, ddds0, dddds0);
        current_x   = x0_nom(:);

        if strcmp(cond_name, 'OnReference')
            mpc.init(current_x, x_ref_pwc, u_ref_pwc);

        elseif strcmp(cond_name, 'Perturbation')
            current_x(1)     = current_x(1) + 0.20;
            current_x(2)     = current_x(2) - 0.20;
            current_x(4:6)   = 0;
            psi0 = s0(4);
            current_x(7)     = cos(psi0 / 2);  % q0
            current_x(8)     = 0;              % q1
            current_x(9)     = 0;              % q2
            current_x(10)    = sin(psi0 / 2);  % q3
            mpc.init(current_x, x_ref_pwc, u_ref_pwc);

        elseif strcmp(cond_name, 'Origin')
            current_x    = zeros(13, 1);
            current_x(3) = -0.02; % 2cm height
            current_x(7) = 1.0;   % q0 = 1 -> no rotation
            if strcmp(scenario_name, 'Lemniscate2')
                current_x(1) = 0.75;
            end
            mpc.init(current_x, x_ref_pwc, u_ref_pwc);
        end

        % Hover beyond reference duration
        x_hover      = x_ref_pwc(:, end);
        x_hover(7)   = 1.0;
        x_hover(4:6) = 0; x_hover(8:10) = 0; x_hover(11:13) = 0;
        u_hover      = [mp.g * mp.m; 0; 0; 0];

        t_steps = 0 : Ts : T_end;
        N_steps = length(t_steps);

        log_states = zeros(N_steps, numel(state_headers));
        log_inputs = zeros(N_steps, numel(input_headers));
        log_errors = zeros(N_steps, numel(error_headers));
        log_perst  = zeros(N_steps, 13);

        sum_sq_pos  = 0;
        sum_sq_full = 0;
        sum_sq_u    = 0;
        solve_times = zeros(N_steps, 1);


        %% RECEDING HORIZON
        for k = 1:N_steps
            t_now = t_steps(k);

            xref_seq = zeros(13, cp.p + 1);
            uref_seq = zeros(4,  cp.p);
            for j = 0:cp.p
                idx = k + j;
                if idx <= size(x_ref_pwc, 2)
                    xref_seq(:, j+1) = x_ref_pwc(:, idx);
                    if j < cp.p, uref_seq(:, j+1) = u_ref_pwc(:, idx); end
                else
                    xref_seq(:, j+1) = x_hover;
                    if j < cp.p, uref_seq(:, j+1) = u_hover; end
                end
            end

            % Receding horizon solver call
            tic;
            [u_opt, ~] = mpc.solve(current_x, xref_seq, uref_seq);
            calc_time  = toc;
            solve_times(k) = calc_time;

            [~, x_sol] = ode45(@(t,x) model.dynamics(t, x, u_opt), ...
                [t_now, t_now + Ts], current_x, options);
            next_x = x_sol(end, :)';

            % Normalize real plant quaternion to prevent integration drift 
            next_x(7:10) = next_x(7:10) / norm(next_x(7:10));

            % Error metrics evaluation
            ref_now  = xref_seq(:, 1);
            uref_now = uref_seq(:, 1);
            e = current_x - ref_now;

            err_pos  = norm(e(1:3));
            err_vel  = norm(e(4:6));
            
            % Exact attitude and yaw extraction for error calculations via Quat_to_euler
            [phi_curr, theta_curr, psi_curr] = Quat_to_euler(current_x(7:10));
            [phi_ref,  theta_ref,  psi_ref ] = Quat_to_euler(ref_now(7:10));
            
            e_att_rad = [phi_curr - phi_ref; theta_curr - theta_ref; psi_curr - psi_ref];
            e_att_rad(3) = atan2(sin(e_att_rad(3)), cos(e_att_rad(3))); % wrap yaw error
            
            err_att    = norm(e_att_rad) * 180/pi;
            e_yaw_raw  = psi_curr - psi_ref;
            e_yaw_wrap = e_att_rad(3);
            
            err_rate = norm(e(11:13)); 
            err_full = norm(e);
            err_u    = norm(u_opt - uref_now);

            sum_sq_pos  = sum_sq_pos  + err_pos^2;
            sum_sq_full = sum_sq_full + err_full^2;
            sum_sq_u    = sum_sq_u    + err_u^2;

            in_track = double(err_pos < pos_threshold);

            log_states(k,:) = [k, current_x', ref_now'];
            log_inputs(k,:) = [k, u_opt', uref_now'];
            log_errors(k,:) = [k, err_pos, err_att, err_vel, err_rate, ...
                e_yaw_raw, e_yaw_wrap, err_full, ...
                sum_sq_pos/k, sqrt(sum_sq_pos/k), ...
                sum_sq_full/k, sqrt(sum_sq_full/k), ...
                sum_sq_u/k, sqrt(sum_sq_u/k), ...
                in_track, calc_time];

            e_ps = abs(e); % raw absolute error per state
            log_perst(k,:) = e_ps';

            current_x = next_x;
        end


        %% POST-LOOP ANALYSIS
        in_track_vec = logical(log_errors(:, 15));

        k_conv = NaN;
        for ki = 1:N_steps-4
            if all(in_track_vec(ki:ki+4)), k_conv = ki; break; end
        end
        t_conv = NaN;
        if ~isnan(k_conv), t_conv = t_steps(k_conv); end

        if sum(in_track_vec) > 10
            rmse_ss_pos  = sqrt(mean(log_errors(in_track_vec, 2).^2));
            rmse_ss_full = sqrt(mean(log_errors(in_track_vec, 8).^2));
        else
            rmse_ss_pos = NaN; rmse_ss_full = NaN;
        end

        max_pos_err = max(log_errors(:,2));
        max_yaw_err = max(abs(log_errors(:,6)));
        max_att_err = max(log_errors(:,3));
        avg_solve   = mean(solve_times) * 1e3;
        worst_solve = max(solve_times)  * 1e3;
        std_solve   = std(solve_times)  * 1e3;
        n_slow      = sum(solve_times > Ts * 0.8);

        perstate_rmse = zeros(1,13);
        for si = 1:13, perstate_rmse(si) = sqrt(mean(log_perst(:,si).^2)); end

        fprintf('  steps=%d  avg=%.1fms  worst=%.1fms  slow=%d\n', ...
            N_steps, avg_solve, worst_solve, n_slow);
        fprintf('  RMSE_pos=%.4fm  RMSE_ss=%.4fm  t_conv=%.2fs\n\n', ...
            log_errors(end,10), rmse_ss_pos, t_conv);


        %% FIGURES
        fig_title = sprintf('%s — %s', scenario_name, cond_name);

        %% Figure 1: State Tracking (13 subplots in 5x3 grid)
        fig1 = figure('Name', sprintf('States: %s', fig_title), ...
            'Color','w', 'Position',[50 50 1200 950]);
        for i = 1:13
            subplot(5,3,i);
            sc = 1; 
            plot(t_steps, log_states(:,i+1)*sc,  'b',   'LineWidth',1.5); hold on; grid on;
            plot(t_steps, log_states(:,i+14)*sc, 'r--', 'LineWidth',1.2);
            ylabel(state_labels{i}, 'Interpreter','latex');
            xlim([0 T_end]);
            if i >= 11, xlabel('t [s]'); end
            if i == 1, legend('Simulated','Reference','Location','best','Interpreter','none'); end
        end
        sgtitle(sprintf('States Tracking: %s', fig_title), 'Interpreter','none');
        saveas(fig1, fullfile(run_dir, 'plot_states.png')); 


        %% Figure 2: Applied Inputs
        fig2 = figure('Name', sprintf('Inputs: %s', fig_title), ...
            'Color','w', 'Position',[150 100 800 600]);
        u_ylims = [cp.u_min(1), cp.u_max(1);
                   cp.u_min(2), cp.u_max(2);
                   cp.u_min(3), cp.u_max(3);
                   cp.u_min(4), cp.u_max(4)];
        for i = 1:4
            subplot(2,2,i);
            plot(t_steps, log_inputs(:,i+1), 'b',   'LineWidth',1.5); hold on; grid on;
            plot(t_steps, log_inputs(:,i+5), 'r--', 'LineWidth',1.2);
            yline(u_ylims(i,1), 'k:', 'LineWidth',0.8);
            yline(u_ylims(i,2), 'k:', 'LineWidth',0.8);
            ylim(u_ylims(i,:));
            ylabel(input_labels{i}); xlim([0 T_end]);
            if i >= 3, xlabel('t [s]'); end
            if i == 1, legend('Applied','Reference','Location','best','Interpreter','none'); end
        end
        sgtitle(sprintf('Inputs Applied: %s', fig_title), 'Interpreter','none');
        saveas(fig2, fullfile(run_dir, 'plot_inputs.png')); 


        %% Figure 3: Tracking Analysis
        fig3 = figure('Name', sprintf('Analysis: %s', fig_title), ...
            'Color','w', 'Position',[250 100 1000 800]);

        subplot(3,2,1);
        semilogy(t_steps, log_errors(:,2)+1e-12, 'b', 'LineWidth',1.5); hold on; grid on;
        yline(pos_threshold, 'r--', sprintf('%.0f cm', pos_threshold*100), ...
            'LabelVerticalAlignment','bottom');
        if ~isnan(t_conv)
            xline(t_conv, 'k-', sprintf('t_{conv} = %.1fs', t_conv), ...
                'LabelOrientation','horizontal');
        end
        ylabel('Position error [m]'); xlabel('t [s]');
        title('Position Error');

        subplot(3,2,2);
        plot(t_steps, log_errors(:,3), 'Color',[0.8 0.2 0.2], 'LineWidth',1.5); grid on;
        ylabel('Attitude error [deg]'); xlabel('t [s]');
        title('Attitude Error (\phi, \theta, \psi)');

        subplot(3,2,3);
        plot(t_steps, log_errors(:,6)*180/pi, 'Color',[0.6 0.6 0.6], 'LineWidth',1.0); hold on; grid on;
        plot(t_steps, log_errors(:,7)*180/pi, 'Color',[0.1 0.4 0.8], 'LineWidth',1.8);
        legend('algebraic','corrected','Location','best');
        ylabel('[deg]'); xlabel('t [s]');
        title('Yaw Error');

        subplot(3,2,4);
        plot(t_steps, log_errors(:,10), 'b',   'LineWidth',1.5, 'DisplayName','Position'); hold on; grid on;
        plot(t_steps, log_errors(:,12), 'm--', 'LineWidth',1.5, 'DisplayName','Full state (13)');
        legend('Location','best');
        ylabel('RMSE'); xlabel('t [s]');
        title('Cumulative RMSE');

        subplot(3,2,5);
        c_pos  = [0.25 0.45 0.80];
        c_vel  = [0.20 0.65 0.30];
        c_quat = [0.85 0.50 0.10];
        c_rate = [0.75 0.20 0.20];
        hold on; grid on;
        h1 = bar(1:3,   perstate_rmse(1:3),   'FaceColor', c_pos);
        h2 = bar(4:6,   perstate_rmse(4:6),   'FaceColor', c_vel);
        h3 = bar(7:10,  perstate_rmse(7:10),  'FaceColor', c_quat);
        h4 = bar(11:13, perstate_rmse(11:13), 'FaceColor', c_rate);
        xticks(1:13);
        ax5 = gca;
        ax5.TickLabelInterpreter = 'latex';
        ax5.XTickLabel = { ...
            '$x$', '$y$', '$z$', ...
            '$\dot{x}_B$', '$\dot{y}_B$', '$\dot{z}_B$', ...
            '$q_0$', '$q_1$', '$q_2$', '$q_3$', ...
            '$p$', '$q_r$', '$r$'};
        xtickangle(30);
        ylabel('RMSE'); xlabel('State');
        title('Per-State RMSE');
        legend([h1 h2 h3 h4], {'position','velocity','quaternions','rates'}, ...
            'Location','northeast','FontSize',7);

        subplot(3,2,6);
        plot(t_steps, solve_times*1e3, 'Color',[0.3 0.3 0.3], 'LineWidth',0.8); hold on; grid on;
        yline(Ts*1e3,     'r--', sprintf('T_s = %.0f ms', Ts*1e3), ...
            'LabelHorizontalAlignment','left');
        yline(Ts*0.8*1e3, 'r:',  '80% budget', ...
            'LabelHorizontalAlignment','right');
        ylabel('Solve time [ms]'); xlabel('t [s]');
        title(sprintf('QP Solve Time  (worst = %.1f ms)', worst_solve));

        sgtitle(sprintf('Tracking Analysis: %s', fig_title), 'Interpreter','none');
        saveas(fig3, fullfile(run_dir, 'plot_analysis.png')); 


        %% GIF - Generate Euler equivalent representation for 3D animator 
        log_states_euler = zeros(N_steps, 25);
        log_states_euler(:, 1) = log_states(:, 1);
        log_states_euler(:, 2:7) = log_states(:, 2:7);
        log_states_euler(:, 11:13) = log_states(:, 12:14);
        log_states_euler(:, 14:19) = log_states(:, 15:20);
        log_states_euler(:, 23:25) = log_states(:, 25:27);
        
        for ki = 1:N_steps
            [phi_s, theta_s, psi_s] = Quat_to_euler(log_states(ki, 8:11)');
            log_states_euler(ki, 8:10) = [phi_s, theta_s, psi_s];
            
            [phi_r, theta_r, psi_r] = Quat_to_euler(log_states(ki, 21:24)');
            log_states_euler(ki, 20:22) = [phi_r, theta_r, psi_r];
        end
        
        N_pwc = size(x_ref_pwc, 2);
        x_ref_pwc_euler = zeros(12, N_pwc);
        x_ref_pwc_euler(1:6, :) = x_ref_pwc(1:6, :);
        x_ref_pwc_euler(10:12, :) = x_ref_pwc(11:13, :);
        for ki = 1:N_pwc
            [phi_r, theta_r, psi_r] = Quat_to_euler(x_ref_pwc(7:10, ki));
            x_ref_pwc_euler(7:9, ki) = [phi_r; theta_r; psi_r];
        end

        gif_name  = fullfile(run_dir, 'tracking.gif');
        GenerateTrackingSimulationGIF(log_states_euler, x_ref_pwc_euler, gif_name, fig_title, '3D', mp);


        %% CSV EXPORTS
        writetable(array2table(log_states, 'VariableNames', state_headers), ...
            fullfile(run_dir, 'states.csv'));
        writetable(array2table(log_inputs, 'VariableNames', input_headers), ...
            fullfile(run_dir, 'inputs.csv'));
        writetable(array2table(log_errors, 'VariableNames', error_headers), ...
            fullfile(run_dir, 'errors.csv'));
        ps_names = arrayfun(@(i) sprintf('state_%d_error', i), 1:13, 'UniformOutput', false);
        writetable(array2table(log_perst, 'VariableNames', ps_names), ...
            fullfile(run_dir, 'per_state_errors.csv'));


        %% ACCUMULATE METRICS
        metrics_summary(end+1,:) = {cond_name, scenario_name, ...
            log_errors(end,10), log_errors(end,12), log_errors(end,14), ...
            avg_solve, N_steps, T_end};

        metrics_extended(end+1,:) = {cond_name, scenario_name, ...
            log_errors(end,10), rmse_ss_pos, ...
            log_errors(end,12), rmse_ss_full, ...
            log_errors(end,9),  log_errors(end,11), ...
            t_conv, ...
            max_pos_err, max_yaw_err*180/pi, max_att_err, ...
            log_errors(end,14), log_errors(end,13), ...
            avg_solve, worst_solve, std_solve, n_slow, ...
            N_steps, T_end};

    end
end


%% EXPORT SUMMARY
cols_sum = {'initial_condition','scenario', ...
    'RMSE_position_m','RMSE_fullstate','RMSE_input', ...
    'avg_solve_time_ms','steps','simulation_duration_s'};
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
    'steps_exceeding_80pct_budget','total_steps','simulation_duration_s'};
writetable(cell2table(metrics_extended, 'VariableNames', cols_ext), ...
    fullfile(fullRunDir, 'MPC_metrics_extended.csv'));

fprintf('Done. Results: %s\n', fullRunDir);

