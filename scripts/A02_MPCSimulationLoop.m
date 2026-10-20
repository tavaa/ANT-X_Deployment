% A02_MPCSimulationLoop.m
% LTV-MPC closed-loop simulation, A02 model (12 states, inertial NED
% velocities, u = [T; L; M; N]), MATLAB environment.
%
%
% Models:  M1 = A02 (flatness / NLP reference) (REFERENCE GENERATOR)
%          M2 = A02 plant, ode45 (replaced by the real drone in the lab)
%          (PLANT)
%          M3 = A02 inside the MPC (prediction, linearization) (PREDICTION
%          MODEL)
%
% Shapes: Circle | Lemniscate | Spiral | Lissajous   
% IC:     OnReference | Perturbation | Origin         
%
% PIPELINE of tick k, between t_k = (k-1)*Ts and t_{k+1} = k*Ts
%   1) APPLY    the command committed at the previous tick (hover at k = 1);
%               it acts on [t_k, t_{k+1}).
%   2) MEASURE  x(t_k): perfect and instantaneous (lab: DroneCmd getters).
%   3) PREDICT  x_hat(t_{k+1}) = Phi(x_k, u_k, Ts) with the A02 model (ode45).
%   4) MPC      linearize along the nominal trajectory from x_hat(t_{k+1}),
%               ZOH discretization, build the QP, solve (OSQP). The returned
%               u* is the command for [t_{k+1}, t_{k+2}).
%   5) DEADLINE the command must reach the drone BEFORE t_{k+1}. If
%               calc time + transmission latency < Ts and the QP is solved,
%               u* is committed and applied at t_{k+1}. Otherwise u* is
%               discarded and the previous command is held one more step.
%   Reference window for the MPC starts at node k+1.
%   The plant evolves over [t_k, t_{k+1}] with fixed input.
%
% Initial input: u_hover = [m*g; 0; 0; 0].
% IC:
%   OnReference  x = x_ref(0); first command = u_ref(0)    (simulation only)
%   Perturbation hover at x_ref(0) + [0.20; -0.20; 0] m, same z, v = 0,
%                angles = 0, rates = 0; first command = u_hover
%   Origin       hover at [0; 0; -1] m, all the rest 0; first command = u_hover
% After the reference (Perturbation, Origin): hover tail at the last point.

clear;
clc;
close all;

%% PATHS
addpath(fullfile('models'));
addpath(fullfile('flatness'));
addpath(fullfile('configs'));
addpath(fullfile('control/MPC'));
addpath(fullfile('trajectories'));
addpath(fullfile('utils'));
addpath(fullfile('casadi'));


%% RUN SELECTION
% Choose Shape to Run
run_shapes = {'Circle'};         % {'Circle','Lemniscate','Spiral','Lissajous'}
run_ICs    = {'Origin'};   % {'OnReference','Perturbation','Origin'}

%% PATHS
[scriptDir, ~, ~] = fileparts(mfilename('fullpath'));
projectRoot = fileparts(scriptDir);

%% OUTPUT DIRECTORY
resultsBaseDir = fullfile(projectRoot, 'results', 'MPC_A02');

if ~exist(resultsBaseDir, 'dir')
    mkdir(resultsBaseDir);
end

todayStr = datestr(now, 'yyyy-mm-dd');
runID = 1;

while true
    runDirName = sprintf('%s_Run_%02d', todayStr, runID);
    fullRunDir = fullfile(resultsBaseDir, runDirName);

    if ~exist(fullRunDir, 'dir')
        mkdir(fullRunDir);
        break;
    end

    runID = runID + 1;
end

fprintf('MPC Simulation Loop, ANT-X (A02, 12-state)\nOutput: %s\n\n', runDirName);


%% SETTINGS
pos_threshold        = 0.05;    % [m] convergence detection threshold
hover_tail_duration  = 3.0;     % [s] hover tail after the reference (Perturbation, Origin)
transmission_latency = 0.01;   % [s] GCS -> drone command latency (0 = ideal) % remove in deployment
make_gif             = true;

% setpoint before IC
pert_xy    = [0.20; -0.20];     % [m] Perturbation offset w.r.t. x_ref(0) (r_n, r_e)
origin_pos = [0; 0; -1];        % [m] Origin start (r_n, r_e, r_d), all shapes

REFERENCE_OFFSET = 1;           % reference window starts at node k+1 (one-step actuation)


%% SYSTEM
mp = ModelParameters();
cp = ControlParameters();
tp = TrajectoryParameters();
Ts = cp.Ts;

plant_ode = odeset('RelTol', 1e-10, 'AbsTol', 1e-12);

plant = A02_SimplifiedModel_Deployment(mp);   % M2: plant (LAB: the real drone)
model = A02_SimplifiedModel_Deployment(mp);   % M3: internal model of the MPC
mpc   = A02_MPCController(model);

Ainv    = mpc.Ainv;                           % [T;L;M;N] -> [T1..T4]
u_hover = [mp.m * mp.g; 0; 0; 0];             % [T; L; M; N]
T_lo    = ControlParameters.u_min_16;         % per-rotor limits [N] (same as the NLP)
T_hi    = ControlParameters.u_max_16;
T_hover_pm = mp.m * mp.g / 4;                 % [N] hover thrust per rotor


%% OSQP warm-up (throwaway hover problem, result discarded)
x_warmup    = zeros(12, 1);
xref_warmup = zeros(12, cp.p + 2);
uref_warmup = repmat(u_hover, 1, cp.p + 2);

mpc.init(x_warmup, xref_warmup, uref_warmup);
mpc.solve(x_warmup, xref_warmup(:, 1:cp.p+1), uref_warmup(:, 1:cp.p));


%% LOG HEADERS
state_headers = {'step', 'x_I','y_I','z_I','v_x','v_y','v_z','phi','theta','psi','p','q','r', ...
                 'ref_xI','ref_yI','ref_zI','ref_vx','ref_vy','ref_vz','ref_phi','ref_theta','ref_psi','ref_p','ref_q','ref_r'};

input_headers = {'step', 'T','L','M','N', 'T_ref','L_ref','M_ref','N_ref', 'T1','T2','T3','T4', 'T1_ref','T2_ref','T3_ref','T4_ref'};

error_headers = {'step', 'position_error_m','attitude_error_deg','velocity_error_ms','rate_error_rads', ...
                 'yaw_error_raw_rad','yaw_error_corrected_rad','fullstate_error', ...
                 'cumulative_mse_pos','cumulative_rmse_pos','cumulative_mse_full','cumulative_rmse_full', ...
                 'cumulative_mse_input','cumulative_rmse_input', 'in_track', ...
                 'solve_time_s','deadline_miss','tx_time_s','total_latency_s','qp_fail'};

state_labels = {'$x_I$ [m]','$y_I$ [m]','$z_I$ [m]', '$v_x$ [m/s]','$v_y$ [m/s]','$v_z$ [m/s]', ...
                '$\phi$ [deg]','$\theta$ [deg]','$\psi$ [deg]', '$p$ [rad/s]','$q$ [rad/s]','$r$ [rad/s]'};

input_labels = {'$T$ [N]','$L$ [N m]','$M$ [N m]','$N$ [N m]'};

angle_idx = [7 8 9];

metrics_summary  = {};
metrics_extended = {};
omega_runs       = {};


%% MAIN: SHAPE -> IC
for s = 1:numel(run_shapes)
    shape = run_shapes{s};

    % Offline PWC reference + validation
    [traj_obj, T_end_base, periodic] = makeShape(shape, tp);
    
    % Generate PWC Trajectory
    nlp_opts = struct('periodic', periodic, 'print_level', 0);
    [x_ref_pwc, u_ref_pwc, t_ref, info_ref] = A02_GeneratePWC_reference(traj_obj, T_end_base, Ts, nlp_opts);
    
    % Validate Reference Trajectory
    rep = A02_Validate_Trajectory(x_ref_pwc, u_ref_pwc, t_ref, info_ref, shape, false);

    if ~rep.pass
        warning('A02_MPCSimulationLoop:RefFail', '%s: reference validation FAILED (simulation continues).', shape);
    end

    % Hover target at the end of the reference
    x_hover_pt = x_ref_pwc(:, end);
    x_hover_pt(4:6)   = 0;
    x_hover_pt(7:8)   = 0;
    x_hover_pt(10:12) = 0;

    for c = 1:numel(run_ICs)
        cond_name = run_ICs{c};

        run_dir = fullfile(fullRunDir, shape, cond_name);

        if ~exist(run_dir, 'dir')
            mkdir(run_dir);
        end

        fig_title = sprintf('%s, %s', shape, cond_name);
        fprintf('%s | %s\n', shape, cond_name);

        % Reference actually used (hover tail for Perturbation / Origin)
        % This can be removed in the actual deployment
        if strcmp(cond_name, 'OnReference')
            x_ref_use = x_ref_pwc;
            u_ref_use = u_ref_pwc;
            T_end_use = T_end_base;
        else % use in the deployment
            n_tail    = round(hover_tail_duration / Ts);
            x_ref_use = [x_ref_pwc, repmat(x_hover_pt, 1, n_tail)];
            u_ref_use = [u_ref_pwc, repmat(u_hover, 1, n_tail)];
            T_end_use = T_end_base + hover_tail_duration;
        end

        % beyond the end of the reference: hover target
        x_hover = x_hover_pt;

        % Initial condition and first command
        current_x = zeros(12, 1);
        
        % Drone setpoints
        % drone.sendPositionSetpoint(IC, psi0);
        switch cond_name
            case 'OnReference' % simulation only -> x_0 = x_ref(t_0)
                current_x = x_ref_pwc(:, 1);
                u_now = u_ref_use(:, 1);
            
            % IC: [x_I(t_0) + 0.20m; y_I(t_0)-0.20m; z_I(t_0); 0, ...,0]
            case 'Perturbation' 
                current_x(1:3) = x_ref_pwc(1:3, 1) + [pert_xy; 0];
                u_now = u_hover;

            % IC: [0; 0; -1.0; 0, ...,0]
            case 'Origin'
                current_x(1:3) = origin_pos;
                u_now = u_hover;

            otherwise
                error('A02_MPCSimulationLoop:IC', 'Unknown IC "%s".', cond_name);
        end
        
        % init controller
        mpc.init(current_x, x_ref_use(:, 1+REFERENCE_OFFSET:end), u_ref_use(:, 1+REFERENCE_OFFSET:end));

        % Grids and logs
        N_intervals      = round(T_end_use / Ts);
        N_steps          = N_intervals + 1;
        N_intervals_base = round(T_end_base / Ts);
        N_steps_base     = N_intervals_base + 1;
        t_steps          = (0:N_steps-1) * Ts;

        log_states = zeros(N_steps, numel(state_headers));
        log_inputs = zeros(N_steps, numel(input_headers));
        log_errors = NaN(N_steps, numel(error_headers));
        log_perst  = zeros(N_steps, 12);

        sum_sq_pos  = 0;
        sum_sq_full = 0;
        sum_sq_u    = 0;

        solve_times     = NaN(N_intervals, 1);
        deadline_miss   = false(N_intervals, 1);
        qp_fail_log     = false(N_intervals, 1);
        tx_times        = NaN(N_intervals, 1);
        total_latencies = NaN(N_intervals, 1);


        % RECEDING HORIZON LOOP - iterate k intervals for the tracking
        for k = 1:N_intervals
            t_k    = t_steps(k);
            t_next = t_k + Ts;

            % APPLY: command committed at the previous tick, acts on [t_k, t_{k+1})
            % LAB: send it to the drone here (drone.sendLMNFsetpointOL());
            % Normalize before send
            u_applied = u_now;

            % MEASURE STATE: perfect and instantaneous
            % LAB: DroneCmd getters
            [x_meas, t_meas] = measureState(current_x, t_k);

            % PREDICT NEXT STATE + MPC Solve
            % (window of the reference starts at node k+1)

            % Obtain reference window
            [xref_seq, uref_seq] = refWindow(x_ref_use, u_ref_use, k + REFERENCE_OFFSET, cp.p, x_hover, u_hover);
            
            % remaining time to compute input = (t_next - t_meas)
            seg = [t_next - t_meas; u_applied];   % u_applied acts from the measurement to t_{k+1}
            
            % call the solver
            tic;
            [u_opt, dbg] = mpc.solve(x_meas, xref_seq, uref_seq, seg);
            calc_time = toc;
            
            % log time to solve
            solve_times(k) = calc_time;

            qp_failed = ~strcmp(dbg.status, 'solved');
            qp_fail_log(k) = qp_failed;

            % DEADLINE: the command must reach the drone before t_{k+1}
            tx_time = transmission_latency;
            total_latency = (t_meas - t_k) + calc_time + tx_time;   % from t_k to arrival on the drone
            tx_times(k) = tx_time;
            total_latencies(k) = total_latency;
            
            % check deadline missed
            deadline_missed = total_latency >= Ts;
            deadline_miss(k) = deadline_missed;

            if deadline_missed
                warning('MPC:ActuationOverrun', 'calc + tx = %.1f ms >= Ts = %.1f ms at t = %.2f s: u_opt discarded, previous command held.', total_latency*1e3, Ts*1e3, t_k);
            end
            
            % handle deadline missed
            if deadline_missed || qp_failed
                u_now = u_applied;   % hold
            else
                u_now = u_opt;       % applied at t_{k+1}
            end

            % PLANT over [t_k, t_{k+1}] with the applied command (ZOH)
            % LAB: the drone evolves by itself; wait until t_{k+1}
            [~, x_ode] = ode45(@(t, x) plant.dynamics(t, x, u_applied), [t_k, t_next], current_x, plant_ode);
            next_x = x_ode(end, :)';

            % Logging (reference at the CURRENT instant t_k)
            [ref_now, uref_now] = refAt(x_ref_use, u_ref_use, k, x_hover, u_hover);

            [e, err] = trackingErrors(current_x, ref_now);
            err_u = norm(Ainv * (u_applied - uref_now));   % input error in rotor-thrust space [N]

            sum_sq_pos  = sum_sq_pos  + err.pos^2;
            sum_sq_full = sum_sq_full + err.full^2;
            sum_sq_u    = sum_sq_u    + err_u^2;

            in_track = double(err.pos < pos_threshold);

            log_states(k, :) = [k, current_x', ref_now'];
            log_inputs(k, :) = [k, u_applied', uref_now', (Ainv*u_applied)', (Ainv*uref_now)'];

            log_errors(k, :) = [k, err.pos, err.att, err.vel, err.rate, err.yaw_raw, err.yaw_wrap, err.full, ...
                                sum_sq_pos/k, sqrt(sum_sq_pos/k), sum_sq_full/k, sqrt(sum_sq_full/k), ...
                                sum_sq_u/k, sqrt(sum_sq_u/k), in_track, ...
                                calc_time, double(deadline_missed), tx_time, total_latency, double(qp_failed)];

            e_ps = abs(e);
            e_ps(angle_idx) = e_ps(angle_idx) * 180/pi;
            log_perst(k, :) = e_ps';

            current_x = next_x;
        end

        % Terminal sample at t = T_end_use (no MPC solve there)
        k_final = N_steps;

        [ref_final, uref_final] = refAt(x_ref_use, u_ref_use, k_final, x_hover, u_hover);

        [e_f, err_f] = trackingErrors(current_x, ref_final);

        sum_sq_pos_f  = sum_sq_pos  + err_f.pos^2;
        sum_sq_full_f = sum_sq_full + err_f.full^2;

        log_states(k_final, :) = [k_final, current_x', ref_final'];
        log_inputs(k_final, :) = [k_final, u_applied', uref_final', (Ainv*u_applied)', (Ainv*uref_final)'];

        log_errors(k_final, :) = [k_final, err_f.pos, err_f.att, err_f.vel, err_f.rate, err_f.yaw_raw, err_f.yaw_wrap, err_f.full, ...
                                  sum_sq_pos_f/N_steps, sqrt(sum_sq_pos_f/N_steps), sum_sq_full_f/N_steps, sqrt(sum_sq_full_f/N_steps), ...
                                  sum_sq_u/N_intervals, sqrt(sum_sq_u/N_intervals), double(err_f.pos < pos_threshold), ...
                                  NaN, NaN, NaN, NaN, NaN];

        e_ps = abs(e_f);
        e_ps(angle_idx) = e_ps(angle_idx) * 180/pi;
        log_perst(k_final, :) = e_ps';

        T_states = array2table(log_states, 'VariableNames', state_headers);
        T_inputs = array2table(log_inputs, 'VariableNames', input_headers);
        T_errors = array2table(log_errors, 'VariableNames', error_headers);


        % POST-LOOP ANALYSIS (reference span only)
        in_track_vec = logical(T_errors.in_track(1:N_steps_base));

        k_conv = NaN;

        for ki = 1:N_steps_base-4
            if all(in_track_vec(ki:ki+4))
                k_conv = ki;
                break;
            end
        end

        t_conv = NaN;

        if ~isnan(k_conv)
            t_conv = t_steps(k_conv);
        end

        if sum(in_track_vec) > 10
            rmse_ss_pos  = sqrt(mean(T_errors.position_error_m(in_track_vec).^2));
            rmse_ss_full = sqrt(mean(T_errors.fullstate_error(in_track_vec).^2));
        else
            rmse_ss_pos  = NaN;
            rmse_ss_full = NaN;
        end

        max_pos_err = max(T_errors.position_error_m(1:N_steps_base));
        max_yaw_err = max(abs(T_errors.yaw_error_corrected_rad(1:N_steps_base)));
        max_att_err = max(T_errors.attitude_error_deg(1:N_steps_base));

        solve_ref_times = solve_times(1:min(N_intervals_base, numel(solve_times)));
        solve_ref_times = solve_ref_times(~isnan(solve_ref_times));

        avg_solve   = mean(solve_ref_times) * 1e3;
        worst_solve = max(solve_ref_times) * 1e3;
        std_solve   = std(solve_ref_times) * 1e3;
        n_slow      = sum(solve_ref_times > Ts * 0.8);

        n_deadline_miss    = sum(deadline_miss(1:N_intervals_base));
        deadline_miss_rate = n_deadline_miss / max(1, N_intervals_base);
        n_qp_fail          = sum(qp_fail_log(1:N_intervals_base));

        avg_tx            = mean(tx_times(1:N_intervals_base), 'omitnan') * 1e3;
        avg_total_latency = mean(total_latencies(1:N_intervals_base), 'omitnan') * 1e3;

        perstate_rmse = zeros(1, 12);

        for si = 1:12
            perstate_rmse(si) = sqrt(mean(log_perst(1:N_steps_base, si).^2));
        end

        rmse_pos_ref  = T_errors.cumulative_rmse_pos(N_steps_base);
        rmse_full_ref = T_errors.cumulative_rmse_full(N_steps_base);
        mse_pos_ref   = T_errors.cumulative_mse_pos(N_steps_base);
        mse_full_ref  = T_errors.cumulative_mse_full(N_steps_base);

        input_err_vec  = vecnorm(log_inputs(1:N_intervals_base, 10:13) - log_inputs(1:N_intervals_base, 14:17), 2, 2);   % rotor thrusts [N]
        rmse_input_ref = sqrt(mean(input_err_vec.^2));
        mse_input_ref  = mean(input_err_vec.^2);

        fprintf(['  samples=%d (ref=%d)  intervals=%d  avg=%.1fms  worst=%.1fms  ', 'slow=%d  deadline_miss=%d (%.1f%%)  qp_fail=%d\n'], ...
            N_steps, N_steps_base, N_intervals, avg_solve, worst_solve, n_slow, n_deadline_miss, 100*deadline_miss_rate, n_qp_fail);

        fprintf('  tx=%.1fms  total_latency=%.1fms\n', avg_tx, avg_total_latency);

        fprintf('  RMSE_pos=%.4fm  RMSE_ss=%.4fm  t_conv=%.2fs\n\n', rmse_pos_ref, rmse_ss_pos, t_conv);

        if n_qp_fail > 0
            warning('A02_MPCSimulationLoop:QPFail', '%s | %s: %d unsolved QPs on the reference span (the run may have diverged).', shape, cond_name, n_qp_fail);
        end


        % FIGURE 1: state tracking
        fig1 = figure('Name', sprintf('States: %s', fig_title), 'Color', 'w', 'Position', [50 50 1400 800]);

        for i = 1:12
            subplot(3, 4, i);
            sc = 1;
            if ismember(i, angle_idx)
                sc = 180/pi;
            end
            plot(t_steps, log_states(:, i+1) * sc, 'b', 'LineWidth', 1.5);
            hold on;
            grid on;
            plot(t_steps, log_states(:, i+13) * sc, 'r--', 'LineWidth', 1.2);
            xline(T_end_base, 'k:', 'LineWidth', 1.2);
            if i == 1
                xline(T_end_base, 'k:', 'Tracking end / Hover start', 'LineWidth', 1.2, 'LabelVerticalAlignment', 'middle', 'LabelHorizontalAlignment', 'left', 'FontSize', 8);
            end
            ylabel(state_labels{i}, 'Interpreter', 'latex');
            xlim([0 T_end_use]);
            if i >= 9
                xlabel('t [s]');
            end
            if i == 1
                legend('Simulated', 'Reference', 'Location', 'best', 'Interpreter', 'none');
            end
        end

        sgtitle(sprintf('States Tracking: %s', fig_title), 'Interpreter', 'none');
        saveas(fig1, fullfile(run_dir, 'plot_states.png'));


        % FIGURE 2: applied inputs, full duration. Top: [T L M N]. Bottom: rotor thrusts with limits
        fig2 = figure('Name', sprintf('Inputs: %s', fig_title), 'Color', 'w', 'Position', [100 100 1400 600]);

        for i = 1:4
            subplot(2, 4, i);
            stairs(t_steps, log_inputs(:, 1+i), 'b', 'LineWidth', 1.5);
            hold on;
            grid on;
            stairs(t_steps, log_inputs(:, 5+i), 'r--', 'LineWidth', 1.2);
            xline(T_end_base, 'k:', 'LineWidth', 1.2);
            ylabel(input_labels{i}, 'Interpreter', 'latex');
            xlim([0 T_end_use]);
            if i == 1
                legend('Applied', 'Reference', 'Location', 'best', 'Interpreter', 'none');
            end
            subplot(2, 4, 4+i);
            stairs(t_steps, log_inputs(:, 9+i), 'b', 'LineWidth', 1.5);
            hold on;
            grid on;
            stairs(t_steps, log_inputs(:, 13+i), 'r--', 'LineWidth', 1.2);
            xline(T_end_base, 'k:', 'LineWidth', 1.2);
            yline(T_lo(i), 'k:', 'LineWidth', 0.8);
            yline(T_hi(i), 'k:', 'LineWidth', 0.8);
            yline(T_hover_pm, 'g:', 'LineWidth', 0.8);
            ylim([T_lo(i), T_hi(i)]);
            xlim([0 T_end_use]);
            ylabel(sprintf('$T_%d$ [N]', i), 'Interpreter', 'latex');
            xlabel('t [s]');
        end

        sgtitle(sprintf('Inputs Applied: %s (bottom: rotor thrusts = A^{-1} u)', fig_title), 'Interpreter', 'none');
        saveas(fig2, fullfile(run_dir, 'plot_inputs.png'));


        % FIGURE 3: tracking analysis, reference span only
        t_steps_ref = t_steps(1:N_steps_base);
        solve_ref   = solve_times(1:N_intervals_base);
        lat_ref     = total_latencies(1:N_intervals_base);
        miss_ref    = deadline_miss(1:N_intervals_base);

        fig3 = figure('Name', sprintf('Analysis: %s', fig_title), 'Color', 'w', 'Position', [250 100 1000 900]);
        subplot(3, 2, 1);
        semilogy(t_steps_ref, T_errors.position_error_m(1:N_steps_base) + 1e-12, 'b', 'LineWidth', 1.5);
        hold on;
        grid on;
        yline(pos_threshold, 'r--', sprintf('%.0f cm', pos_threshold*100), 'LabelVerticalAlignment', 'bottom');

        if ~isnan(t_conv) && t_conv <= T_end_base
            xline(t_conv, 'k-', sprintf('t_{conv} = %.1fs', t_conv), 'LabelOrientation', 'horizontal');
        end

        xlim([0 T_end_base]);
        ylabel('Position error [m]');
        xlabel('t [s]');
        title('Position Error (tracking span)');
        subplot(3, 2, 2);
        plot(t_steps_ref, T_errors.attitude_error_deg(1:N_steps_base), 'Color', [0.8 0.2 0.2], 'LineWidth', 1.5);
        grid on;
        xlim([0 T_end_base]);
        ylabel('Attitude error [deg]');
        xlabel('t [s]');
        title('Attitude Error  (\phi, \theta, \psi)');
        subplot(3, 2, 3);
        plot(t_steps_ref, T_errors.yaw_error_raw_rad(1:N_steps_base) * 180/pi, 'Color', [0.6 0.6 0.6], 'LineWidth', 1.0);

        hold on;
        grid on;

        plot(t_steps_ref, T_errors.yaw_error_corrected_rad(1:N_steps_base) * 180/pi, 'Color', [0.1 0.4 0.8], 'LineWidth', 1.8);
        xlim([0 T_end_base]);
        legend('algebraic', 'corrected', 'Location', 'best');

        ylabel('[deg]');
        xlabel('t [s]');
        title('Yaw Error');

        subplot(3, 2, 4);

        plot(t_steps_ref, T_errors.velocity_error_ms(1:N_steps_base), 'Color', [0.2 0.65 0.3], 'LineWidth', 1.5);
        hold on;
        grid on;
        plot(t_steps_ref, T_errors.rate_error_rads(1:N_steps_base), 'Color', [0.75 0.2 0.2], 'LineWidth', 1.5);
        xlim([0 T_end_base]);
        legend('velocity [m/s]', 'body rates [rad/s]', 'Location', 'best');
        xlabel('t [s]');
        title('Velocity and Body-Rate Error');
        subplot(3, 2, 5);

        hold on;
        grid on;

        h1 = bar(1:3,   perstate_rmse(1:3),   'FaceColor', [0.25 0.45 0.80]);
        h2 = bar(4:6,   perstate_rmse(4:6),   'FaceColor', [0.20 0.65 0.30]);
        h3 = bar(7:9,   perstate_rmse(7:9),   'FaceColor', [0.85 0.50 0.10]);
        h4 = bar(10:12, perstate_rmse(10:12), 'FaceColor', [0.75 0.20 0.20]);
        
        xticks(1:12);

        ax5 = gca;
        ax5.TickLabelInterpreter = 'latex';
        ax5.XTickLabel = {'$x_I$', '$y_I$', '$z_I$', '$v_x$', '$v_y$', '$v_z$', '$\phi$', '$\theta$', '$\psi$', '$p$', '$q$', '$r$'};
        xtickangle(30);

        ylabel('RMSE');
        xlabel('State');
        title('Per-State RMSE (tracking span)');

        legend([h1 h2 h3 h4], {'position', 'velocity', 'angles [deg]', 'rates'}, 'Location', 'northeast', 'FontSize', 7);

        subplot(3, 2, 6);

        plot(t_steps_ref(1:end-1), solve_ref * 1e3, 'Color', [0.3 0.3 0.3], 'LineWidth', 0.8);

        hold on;
        grid on;

        plot(t_steps_ref(1:end-1), lat_ref * 1e3, 'Color', [0.1 0.4 0.8], 'LineWidth', 0.8);

        if any(miss_ref)
            plot(t_steps_ref(miss_ref), lat_ref(miss_ref) * 1e3, 'rx', 'MarkerSize', 8);
        end

        yline(Ts * 1e3, 'r--', sprintf('T_s = %.0f ms', Ts*1e3), 'LabelHorizontalAlignment', 'left');

        xlim([0 T_end_base]);

        ylabel('[ms]');
        xlabel('t [s]');
        legend('solve time', 'calc + tx (deadline)', 'Location', 'best');
        title(sprintf('Controller time  (worst = %.1f ms, misses = %d)', worst_solve, n_deadline_miss));

        sgtitle(sprintf('Tracking Analysis: %s (tracking span only, hover excluded)', fig_title), 'Interpreter', 'none');

        saveas(fig3, fullfile(run_dir, 'plot_analysis.png'));


        % GIF
        if make_gif
            try
                GenerateTrackingSimulationGIF(log_states, x_ref_use, fullfile(run_dir, 'tracking.gif'), fig_title, '3D', mp, C00_PlotDrone);
            catch gif_err
                warning('A02_MPCSimulationLoop:GIF', 'GIF not generated: %s', gif_err.message);
            end
        end


        % CSV and MAT, full duration logs
        writetable(T_states, fullfile(run_dir, 'states.csv'));
        writetable(T_inputs, fullfile(run_dir, 'inputs.csv'));
        writetable(T_errors, fullfile(run_dir, 'errors.csv'));

        ps_names = arrayfun(@(i) sprintf('state_%d_error', i), 1:12, 'UniformOutput', false);
        writetable(array2table(log_perst, 'VariableNames', ps_names), fullfile(run_dir, 'per_state_errors.csv'));

        save(fullfile(run_dir, 'run_data.mat'), 'T_states', 'T_inputs', 'T_errors', 'x_ref_use', 'u_ref_use', 't_steps', 'solve_times', 'deadline_miss', 'qp_fail_log');


        % JSON CONFIGURATION
        cfg = struct();
        cfg.run         = struct('shape', shape, 'initial_condition', cond_name, 'date', datestr(now, 'yyyy-mm-dd HH:MM:SS'));
        cfg.model       = struct('plant', class(plant), 'mpc_model', class(model), 'm', mp.m, 'g', mp.g, 'Jxx', mp.Jxx, 'Jyy', mp.Jyy, 'Jzz', mp.Jzz);
        cfg.control     = struct('Ts', Ts, 'horizon', cp.p, 'Q', diag(mpc.Q), 'Qf', diag(mpc.Qf), 'R', diag(mpc.R), ...
                                 'rotor_min', T_lo, 'rotor_max', T_hi, 'x_min', cp.x_min_12, 'x_max', cp.x_max_12, ...
                                 'reference_offset', REFERENCE_OFFSET);
        cfg.timing      = struct('measurement', 'perfect, instantaneous', 'transmission_latency_s', transmission_latency, ...
                                 'deadline', 'calc + tx < Ts, else previous command held');
        cfg.reference   = struct('T_end', T_end_base, 'periodic', periodic, 'nlp', info_ref.opts, 'validation_pass', rep.pass);
        cfg.ic          = struct('pert_xy', pert_xy, 'origin_pos', origin_pos, 'hover_tail_s', hover_tail_duration);
        cfg.trajectory  = traj_obj_params(shape, tp);

        fid = fopen(fullfile(run_dir, 'run_config.json'), 'w');

        if fid ~= -1
            fwrite(fid, jsonencode(cfg, 'PrettyPrint', true), 'char');
            fclose(fid);
        end


        % ACCUMULATE METRICS (reference span only)
        metrics_summary(end+1, :) = {cond_name, shape, rmse_pos_ref, rmse_full_ref, rmse_input_ref, avg_solve, N_steps_base, T_end_base};  %#ok<SAGROW>

        metrics_extended(end+1, :) = {cond_name, shape, rmse_pos_ref, rmse_ss_pos, rmse_full_ref, rmse_ss_full, mse_pos_ref, mse_full_ref, ...
            t_conv, max_pos_err, max_yaw_err*180/pi, max_att_err, rmse_input_ref, mse_input_ref, ...
            avg_solve, worst_solve, std_solve, n_slow, n_deadline_miss, deadline_miss_rate, avg_tx, avg_total_latency, n_qp_fail, ...
            N_steps_base, T_end_base};  %#ok<SAGROW>

        omega_runs{end+1} = struct('name', fig_title, 't', t_steps, 'T', log_inputs(:, 10:13));  %#ok<SAGROW>
    end
end


%% EXPORT SUMMARY
cols_sum = {'initial_condition','scenario', 'RMSE_position_m','RMSE_fullstate','RMSE_input_rotor_N', 'avg_solve_time_ms','steps_reference_span','reference_duration_s'};

writetable(cell2table(metrics_summary, 'VariableNames', cols_sum), fullfile(fullRunDir, 'MPC_metrics_summary.csv'));

cols_ext = {'initial_condition','scenario', 'rmse_position_all_steps_m','rmse_position_steadystate_m', ...
            'rmse_fullstate_all_steps','rmse_fullstate_steadystate', 'mse_position_m2','mse_fullstate', ...
            'convergence_time_s', 'peak_position_error_m','peak_yaw_error_deg','peak_attitude_error_deg', ...
            'rmse_input_rotor_N','mse_input_rotor', 'avg_solve_time_ms','worst_solve_time_ms','std_solve_time_ms', ...
            'steps_exceeding_80pct_budget','deadline_misses','deadline_miss_rate', ...
            'avg_tx_time_ms','avg_total_latency_ms','qp_failures', 'steps_reference_span','reference_duration_s'};

writetable(cell2table(metrics_extended, 'VariableNames', cols_ext), fullfile(fullRunDir, 'MPC_metrics_extended.csv'));

fprintf('Done. Results: %s\n', fullRunDir);


%% DEBUG: rotor angular velocities, Omega_i = sqrt(T_i / K_T)
K_T = mp.K_T;
Omega_min = sqrt(T_lo(1) / K_T);
Omega_max = sqrt(T_hi(1) / K_T);

motor_colors = {[0.20 0.45 0.80], [0.85 0.20 0.20], [0.20 0.65 0.30], [0.60 0.30 0.70]};

fig_om = figure('Name', 'Omega_i plausibility check', 'Color', 'w', 'Position', [50 50 1400 700]);
tl = tiledlayout('flow', 'Padding', 'compact', 'TileSpacing', 'compact');

for r = 1:numel(omega_runs)
    nexttile;

    for i = 1:4
        plot(omega_runs{r}.t, sqrt(max(omega_runs{r}.T(:, i), 0) / K_T), 'Color', motor_colors{i}, 'LineWidth', 1.1);
        hold on;
    end

    grid on;

    yline(Omega_min, 'k:', 'LineWidth', 0.7);
    yline(Omega_max, 'k:', '$\Omega_{\max}$', 'Interpreter', 'latex', 'LineWidth', 0.7, 'LabelHorizontalAlignment', 'right');
    yline(sqrt(T_hover_pm / K_T), 'g:', '$\Omega_{\mathrm{hover}}$', 'Interpreter', 'latex', 'LineWidth', 0.7, 'LabelHorizontalAlignment', 'right');

    title(omega_runs{r}.name, 'FontSize', 9);
    ylabel('$\Omega$ [rad/s]', 'Interpreter', 'latex');
    xlabel('t [s]');
end

lgd = legend({'$\Omega_1$', '$\Omega_2$', '$\Omega_3$', '$\Omega_4$'}, 'Interpreter', 'latex', 'Orientation', 'horizontal');
lgd.Layout.Tile = 'south';

sgtitle(tl, '$\mathrm{Rotor\ angular\ velocities:}\quad \Omega_i = \sqrt{T_i/K_T}$', 'Interpreter', 'latex');

saveas(fig_om, fullfile(fullRunDir, 'plot_omega.png'));


%% LOCAL FUNCTIONS
function [traj_obj, T_end, periodic] = makeShape(name, tp)
    % MAKESHAPE  Trajectory object, duration and periodicity of the shape.
    switch name
        case 'Circle'
            traj_obj = ShapeCircle(tp.Circle);            T_end = tp.T_duration_circle;      periodic = true;
        case 'Lemniscate'
            traj_obj = ShapeLemniscate2(tp.Lemniscate);   T_end = tp.T_duration_lemniscate;  periodic = true;
        case 'Spiral'
            traj_obj = ShapeSpiral(tp.Spiral);            T_end = tp.T_duration_spiral;      periodic = false;   % open curve
        case 'Lissajous'
            traj_obj = ShapeLissajous3D(tp.Lissajous3D);  T_end = tp.T_duration_lissajous3D; periodic = true;
        otherwise
            error('A02_MPCSimulationLoop:Shape', 'Unknown shape "%s".', name);
    end
end


function p = traj_obj_params(name, tp)
    % TRAJ_OBJ_PARAMS  Parameter struct of the shape (for the JSON).
    switch name
        case 'Circle',     p = tp.Circle;
        case 'Lemniscate', p = tp.Lemniscate;
        case 'Spiral',     p = tp.Spiral;
        case 'Lissajous',  p = tp.Lissajous3D;
    end
end


function [x_meas, t_meas] = measureState(x_plant, t_now)
    % MEASURESTATE  Simulated state acquisition: perfect and instantaneous.
    %               The attitude arrives wrapped in [-pi, pi], as from the drone.
    %               LAB: replace with the DroneCmd getters (position, velocity,
    %               Euler angles, body rates) and t_meas = toc(t0).
    %
    %                   [pn, pe, pd] = drone.getPosition();
    %                   [vn, ve, vd] = drone.getVelocity();
    %                   [ph, th, ps] = drone.getEuler();
    %                   [p,  q,  r ] = drone.getAngularVelocity();
    x_meas = x_plant;
    x_meas(9) = atan2(sin(x_meas(9)), cos(x_meas(9)));
    t_meas = t_now;
end


function [e, err] = trackingErrors(x, x_ref)
    % TRACKINGERRORS  Error vector (yaw wrapped) and its norms.
    e = x - x_ref;

    err.yaw_raw  = e(9);
    err.yaw_wrap = atan2(sin(e(9)), cos(e(9)));
    e(9) = err.yaw_wrap;

    err.pos  = norm(e(1:3));
    err.vel  = norm(e(4:6));
    err.att  = norm(e(7:9)) * 180/pi;   % [deg]
    err.rate = norm(e(10:12));
    err.full = norm(e);
end


function [xw, uw] = refWindow(x_ref, u_ref, k0, p, x_hover, u_hover)
    % REFWINDOW  Reference window over the horizon, starting at sample k0.
    %            Beyond the end of the reference: hover target.
    xw = zeros(size(x_ref, 1), p + 1);
    uw = zeros(4, p);
    T = size(x_ref, 2);

    for j = 0:p
        idx = k0 + j;

        if idx <= T
            xw(:, j+1) = x_ref(:, idx);

            if j < p
                uw(:, j+1) = u_ref(:, idx);
            end
        else
            xw(:, j+1) = x_hover;

            if j < p
                uw(:, j+1) = u_hover;
            end
        end
    end
end


function [xr, ur] = refAt(x_ref, u_ref, k, x_hover, u_hover)
    % REFAT  Reference state/input at sample k.
    if k <= size(x_ref, 2)
        xr = x_ref(:, k);
        ur = u_ref(:, k);
    else
        xr = x_hover;
        ur = u_hover;
    end
end