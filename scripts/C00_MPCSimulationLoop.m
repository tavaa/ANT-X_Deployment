% C00_MPCSimulationLoop.m
% LTV-MPC closed-loop simulation for the ANT-X quadrotor, C00 model.
% Scenarios: Circle, Lemniscate, Spiral, Lissajous
% Initial conditions: OnReference, Perturbation, Origin
%
% Pipeline:
%   model      C00_ANTX_quadcopter (16 states, motor lag tau_p)
%   reference  C00_GeneratePWC_reference (NLP, 16 states) + C00_TrajectoryValidation
%   controller C00_MPCController + C00_QPBuilder (LTV-MPC, OSQP)
%
% Timing: fixed one-step actuation is used throughout the simulation.
% The input computed at t_k is intended for [t_{k+1}, t_{k+2}) and the
% reference window starts at k+1. The state the controller receives is
% OLDER than t_k (age a) and noisy: the controller predicts x(t_{k+1}) by
% integrating the model from the measurement time stamp t_k - a with the
% commands that really acted in between, i.e. the HISTORY
%     [a, Ts ; u_prev, u_now]          (u_prev acted on [t_k-a, t_k))
% If calculation + Wi-Fi transmission misses the t_{k+1} deadline, or the QP
% is not solved, the newly computed command is discarded and the previously
% committed command is held for one more step.
%
% Plant = the same C00 model integrated with ode45 (no model mismatch).
% State acquisition can be simulated either as perfect state or as delayed +
% noisy state. In the lab x(1:12) is expected from mocap (position, attitude)
% and the PX4 estimator (inertial velocity, body rates), x(13:16) from ESC
% RPM telemetry (T = K_T*Omega^2).

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

%% OUTPUT DIRECTORY
[scriptDir, ~, ~] = fileparts(mfilename('fullpath'));
projectRoot = fileparts(scriptDir);
resultsBaseDir = fullfile(projectRoot, 'results/MPC_C00/');

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

fprintf('MPC Simulation Loop — ANT-X (C00, 16-state)\nOutput: %s\n\n', runDirName);


%% SETTINGS
pos_threshold = 0.05;  % [m] convergence detection threshold
hover_tail_duration = 3.0;  % [s] hover tail (Perturbation, Origin)
perfect_state = true;  % true: perfect state; false: delayed + noisy state
transmission_latency = 0.010;  % [s] nominal one-way Wi-Fi/GCS command latency
rng_seed = 42;  % reproducible measurement noise
make_gif = true;  % GenerateTrackingSimulationGIF

% [SIM] Age of the data of each measurement source.
% The controller needs ONE time stamp, so every group is aligned to
% the OLDEST source: the age used is the maximum.
state_age_src = struct('mocap', 0.040, 'px4', 0.020, 'esc', 0.030);  % position + attitude | inertial velocity + body rates | thrust (RPM telemetry)

state_acquisition_delay = max([state_age_src.mocap, state_age_src.px4, state_age_src.esc]);

% [SIM] Measurement-noise assumptions for the non-perfect-state case.
state_noise_std = [5e-3; 5e-3; 5e-3; 2e-2; 2e-2; 2e-2; deg2rad(1); deg2rad(1); deg2rad(1); 3e-2; 3e-2; 3e-2; 1e-2; 1e-2; 1e-2; 1e-2];  % position [m] | velocity [m/s] | attitude [rad] | body rates [rad/s] | thrust states [N]

% Origin start [r_n; r_e] per scenario. Lemniscate pass
% through the centre of the room: its Origin start is moved off-centre
origin_xy = struct('Circle', [0.00; 0.00], 'Lemniscate', [0.75; 0.00], 'Spiral', [0.00; 0.00], 'Lissajous', [0.00; 0.00]);


%% SYSTEM
mp = ModelParameters();
cp = ControlParameters();
tp = TrajectoryParameters();
Ts = cp.Ts;
options = odeset('RelTol', 1e-8, 'AbsTol', 1e-9);

assert(state_acquisition_delay <= Ts, 'The history logic ([a, Ts; u_prev, u_now]) assumes a state age not larger than Ts.');

% model and controller
model = C00_ANTX_quadcopter(mp, true);  % delay_motors = true
mpc = C00_MPCController(model);

T_hover_pm = mp.m * mp.g / 4;  % [N] per-motor hover thrust
REFERENCE_OFFSET = 1;  % fixed one-step actuation/reference alignment


%% OSQP warm-up.
% Throwaway generic hover-consistent problem, result discarded.
x_warmup = [zeros(12,1); T_hover_pm*ones(4,1)];
xref_warmup = [zeros(12, cp.p + 2); T_hover_pm*ones(4, cp.p + 2)];
uref_warmup = T_hover_pm * ones(4, cp.p + 2);

mpc.init(x_warmup, xref_warmup, uref_warmup);
mpc.solve(x_warmup, xref_warmup(:, 1:cp.p+1), uref_warmup(:, 1:cp.p));


%% TRAJECTORIES (offline PWC reference generation, NLP, 16-state)
traj_circle = ShapeCircle(tp.Circle);
traj_lemniscate = ShapeLemniscate2(tp.Lemniscate);
traj_spiral = ShapeSpiral(tp.Spiral);
traj_lissajous = ShapeLissajous3D(tp.Lissajous3D);

nlp_closed = struct('periodic', true, 'print_level', 0);
nlp_open = struct('periodic', false, 'print_level', 0);  % spiral: open curve

[x_circle, u_circle, t_circle, i_circle] = C00_GeneratePWC_reference(traj_circle, tp.T_duration_circle, Ts, nlp_closed);

[x_lemniscate, u_lemniscate, t_lemniscate, i_lemniscate] = C00_GeneratePWC_reference(traj_lemniscate, tp.T_duration_lemniscate, Ts, nlp_closed);

[x_spiral, u_spiral, t_spiral, i_spiral] = C00_GeneratePWC_reference(traj_spiral, tp.T_duration_spiral, Ts, nlp_open);

[x_lissajous, u_lissajous, t_lissajous, i_lissajous] = C00_GeneratePWC_reference(traj_lissajous, tp.T_duration_lissajous3D, Ts, nlp_closed);


%% SCENARIOS  {name, traj_obj, T_end, x_ref, u_ref}
scenarios = {
    %'Circle', traj_circle, tp.T_duration_circle, x_circle, u_circle;
    %'Lemniscate', traj_lemniscate, tp.T_duration_lemniscate, x_lemniscate, u_lemniscate;
    'Spiral', traj_spiral, tp.T_duration_spiral, x_spiral, u_spiral;
    'Lissajous', traj_lissajous, tp.T_duration_lissajous3D, x_lissajous, u_lissajous;
};


%% REFERENCE VALIDATION (before closing the loop; console only)
val_in = {x_circle, u_circle, t_circle, i_circle;
    x_lemniscate, u_lemniscate, t_lemniscate, i_lemniscate;
    x_spiral, u_spiral, t_spiral, i_spiral;
    x_lissajous, u_lissajous, t_lissajous, i_lissajous};

for s = 1:size(scenarios, 1)
    rep = C00_TrajectoryValidation(val_in{s,1}, val_in{s,2}, val_in{s,3}, val_in{s,4}, scenarios{s,1}, false);

    if ~rep.pass
        warning('C00_MPCSimulationLoop:RefFail', '%s: reference validation FAILED (simulation continues).', scenarios{s,1});
    end
end


%% JSON CONFIGURATION
cfg.model = struct('name', class(model), 'flatness', 'C00_FlatnessMap', 'reference', 'C00_GeneratePWC_reference (NLP, 16 states)', 'm', mp.m, 'g', mp.g, 'Jxx', mp.Jxx, 'Jyy', mp.Jyy, 'Jzz', mp.Jzz, 'Lp', mp.Lp, 'Mq', mp.Mq, 'Nr', mp.Nr, 'tau_p', mp.tau_p);

cfg.control = struct('Ts', cp.Ts, 'horizon', cp.p, 'Q', cp.Q_diag_16, 'Qf', cp.Qf_diag_16, 'R', cp.R_diag_16, 'u_min', cp.u_min_16, 'u_max', cp.u_max_16, 'x_min', cp.x_min_16, 'x_max', cp.x_max_16, 'reference_offset', REFERENCE_OFFSET);

cfg.timing = struct('state_age_src', state_age_src, 'state_acquisition_delay', state_acquisition_delay, 'transmission_latency', transmission_latency, 'Ts', Ts);

cfg.measurement = struct('perfect_state', perfect_state, 'noise_std', state_noise_std, 'rng_seed', rng_seed);

cfg.trajectories = struct('Circle', tp.Circle, 'Lemniscate', tp.Lemniscate, 'Spiral', tp.Spiral, 'Lissajous', tp.Lissajous3D);

cfg.origin_xy = origin_xy;

fid = fopen(fullfile(fullRunDir, 'simulation_config.json'), 'w');

if fid ~= -1
    fwrite(fid, jsonencode(cfg, 'PrettyPrint', true), 'char');
    fclose(fid);
end


%% LOG HEADERS
state_headers = {'step', 'r_n','r_e','r_d','v_n','v_e','v_d','phi','theta','psi','p','q','r','T1','T2','T3','T4', 'ref_rn','ref_re','ref_rd','ref_vn','ref_ve','ref_vd','ref_phi','ref_theta','ref_psi','ref_p','ref_q','ref_r','ref_T1','ref_T2','ref_T3','ref_T4'};

input_headers = {'step', 'T1','T2','T3','T4', 'T1_ref','T2_ref','T3_ref','T4_ref'};

error_headers = {'step', 'position_error_m','attitude_error_deg','velocity_error_ms','rate_error_rads', 'thrust_error_N', 'yaw_error_raw_rad','yaw_error_corrected_rad','fullstate_error', 'cumulative_mse_pos','cumulative_rmse_pos', 'cumulative_mse_full','cumulative_rmse_full', 'cumulative_mse_input','cumulative_rmse_input', 'in_track','solve_time_s','deadline_miss','tx_time_s','total_latency_s','state_age_s','qp_fail'};

state_labels = {'$r_n$ [m]','$r_e$ [m]','$r_d$ [m]', '$v_n$ [m/s]','$v_e$ [m/s]','$v_d$ [m/s]', '$\phi$ [deg]','$\theta$ [deg]','$\psi$ [deg]', '$p$ [rad/s]','$q$ [rad/s]','$r$ [rad/s]', '$T_1$ [N]','$T_2$ [N]','$T_3$ [N]','$T_4$ [N]'};

input_labels = {'T_1^d [N]','T_2^d [N]','T_3^d [N]','T_4^d [N]'};

angle_idx = [7 8 9];


%% MAIN LOOP
rng(rng_seed);

%initial_conditions = {'OnReference','Perturbation','Origin'};
initial_conditions = {'Perturbation'}; % test single condition

metrics_summary = {};
metrics_extended = {};

for c = 1:numel(initial_conditions)
    cond_name = initial_conditions{c};

    for s = 1:size(scenarios, 1)

        scenario_name = scenarios{s, 1};
        traj_obj = scenarios{s, 2};
        T_end_base = scenarios{s, 3};  % end of the REFERENCE
        x_ref_pwc = scenarios{s, 4};  % [16 x (N+1)] NLP reference, no padding
        u_ref_pwc = scenarios{s, 5};  % [4  x (N+1)] last column = padding

        run_dir = fullfile(fullRunDir, scenario_name, cond_name);

        if ~exist(run_dir, 'dir')
            mkdir(run_dir);
        end

        fprintf('%s | %s\n', scenario_name, cond_name);

        %% Hover at the end of the reference (Perturbation / Origin)
        x_hover_pt = x_ref_pwc(:, end);

        x_hover_pt(4:6) = 0;
        x_hover_pt(7:8) = 0;
        x_hover_pt(10:12) = 0;
        x_hover_pt(13:16) = T_hover_pm;

        if strcmp(cond_name, 'Perturbation') || strcmp(cond_name, 'Origin')
            n_tail = round(hover_tail_duration / Ts);
            x_ref_use = [x_ref_pwc, repmat(x_hover_pt, 1, n_tail)];
            u_ref_use = [u_ref_pwc, repmat(T_hover_pm*ones(4,1), 1, n_tail)];
            T_end_use = T_end_base + hover_tail_duration;
        else
            x_ref_use = x_ref_pwc;
            u_ref_use = u_ref_pwc;
            T_end_use = T_end_base;
        end

        % beyond the end of the reference: hover target
        x_hover = x_hover_pt;
        u_hover = T_hover_pm * ones(4,1);

        %% Initial condition
        current_x = x_ref_pwc(:, 1);  % = x^ref(0), fixed by the NLP

        if strcmp(cond_name, 'Perturbation')
            current_x(1) = current_x(1) + 0.20;
            current_x(2) = current_x(2) - 0.20;
            current_x(4:6) = 0;
            current_x(7:8) = 0;
            current_x(9) = 0;
            current_x(10:12) = 0;

        elseif strcmp(cond_name, 'Origin')
            current_x = zeros(16, 1);
            current_x(1:2) = origin_xy.(scenario_name);
            current_x(3) = -0.02;
            current_x(13:16) = cp.u_min_16;
        end

        % input acting on the first interval: reference input if the drone
        % starts on the reference, hover thrust otherwise
        if strcmp(cond_name, 'OnReference')
            u_now = u_ref_use(:, 1);
        else
            u_now = T_hover_pm * ones(4,1);
        end

        u_prev = u_now;  % command that acted before t_1 (history for the first ticks)

        mpc.init(current_x, x_ref_use(:, 1+REFERENCE_OFFSET:end), u_ref_use(:, 1+REFERENCE_OFFSET:end));

        % N_steps counts samples; N_intervals counts plant propagation intervals.
        N_intervals = round(T_end_use / Ts);
        N_steps = N_intervals + 1;
        N_intervals_base = round(T_end_base / Ts);
        N_steps_base = N_intervals_base + 1;
        t_steps = (0:N_steps-1) * Ts;

        log_states = zeros(N_steps, numel(state_headers));
        log_inputs = zeros(N_steps, numel(input_headers));
        log_errors = NaN(N_steps, numel(error_headers));
        log_perst = zeros(N_steps, 16);

        % DENSE history of the TRUE plant state (ode45 output), needed to
        % emulate a measurement taken at an older time stamp.
        hist_t = 0;
        hist_x = current_x;

        sum_sq_pos = 0;
        sum_sq_full = 0;
        sum_sq_u = 0;

        solve_times = NaN(N_intervals, 1);
        deadline_miss = false(N_intervals, 1);
        qp_fail_log = false(N_intervals, 1);
        tx_times = NaN(N_intervals, 1);
        total_latencies = NaN(N_intervals, 1);


        %% RECEDING HORIZON
        for k = 1:N_intervals
            t_now = t_steps(k);

            % State acquisition model. In the realistic simulation case the
            % controller receives a stale/noisy state, not current_x.
            if perfect_state
                x_meas = current_x;
                state_age = 0;
                seg = [Ts; u_now];  % history: only the committed command
            else
                t_meas = max(0, t_now - state_acquisition_delay);
                x_meas = delayedState(hist_t, hist_x, t_meas);
                x_meas = x_meas + state_noise_std .* randn(16,1);
                x_meas(9) = atan2(sin(x_meas(9)), cos(x_meas(9)));  % attitude arrives wrapped in [-pi,pi]
                state_age = t_now - t_meas;  % actual age (shorter in the first ticks)
                seg = [state_age, Ts; u_prev, u_now];  % history: u_prev on [t_k-a,t_k), u_now on [t_k,t_k+Ts)
            end

            % One-step actuation is mandatory: u_opt is intended for the NEXT interval.
            [xref_seq, uref_seq] = refWindow(x_ref_use, u_ref_use, k + REFERENCE_OFFSET, cp.p, x_hover, u_hover);

            tic;
            [u_opt, dbg] = mpc.solve(x_meas, xref_seq, uref_seq, seg);
            calc_time = toc;

            solve_times(k) = calc_time;

            qp_failed = ~strcmp(dbg.status, 'solved');  % OSQP not solved: do not trust u_opt
            qp_fail_log(k) = qp_failed;

            % [SIM] nominal one-way Wi-Fi/GCS command transmission latency.
            % The acquisition delay happened before this controller execution,
            % therefore it is not added to the command deadline.
            tx_time = transmission_latency;
            total_actuation_delay = calc_time + tx_time;
            tx_times(k) = tx_time;
            total_latencies(k) = total_actuation_delay;

            deadline_missed = total_actuation_delay >= Ts;
            deadline_miss(k) = deadline_missed;

            if deadline_missed
                warning('MPC:ActuationOverrun', 'calc+Wi-Fi latency = %.1f ms >= Ts = %.1f ms at t = %.2f s: u_opt discarded; previous command is held.', total_actuation_delay*1e3, Ts*1e3, t_now);
            end

            % A missed deadline or an unsolved QP leaves the committed command unchanged.
            if deadline_missed || qp_failed
                u_next_committed = u_now;
            else
                u_next_committed = u_opt;
            end

            % The command already committed before t_k acts over the whole
            % current interval [t_k, t_{k+1}).
            [t_ode, x_ode] = ode45(@(t,x) model.dynamics(t, x, u_now), [t_now, t_now + Ts], current_x, options);

            next_x = x_ode(end, :)';

            hist_t = [hist_t, t_ode(2:end).'];        %#ok<AGROW>
            hist_x = [hist_x, x_ode(2:end, :).'];     %#ok<AGROW>

            u_applied = u_now;  % input actually applied on [t_k,t_{k+1})
            u_prev = u_applied;  % it is the "previous" command at the next tick
            u_now = u_next_committed;  % command committed for [t_{k+1},t_{k+2})

            %% Logging (reference at the CURRENT instant t_k)
            [ref_now, uref_now] = refAt(x_ref_use, u_ref_use, k, x_hover, u_hover);

            e = current_x - ref_now;
            e_yaw_raw = e(9);
            e_yaw_wrap = atan2(sin(e_yaw_raw), cos(e_yaw_raw));
            e(9) = e_yaw_wrap;

            err_pos = norm(e(1:3));
            err_vel = norm(e(4:6));
            err_att = norm(e(7:9)) * 180/pi;
            err_rate = norm(e(10:12));
            err_thrust = norm(e(13:16));
            err_full = norm(e);
            err_u = norm(u_applied - uref_now);

            sum_sq_pos = sum_sq_pos + err_pos^2;
            sum_sq_full = sum_sq_full + err_full^2;
            sum_sq_u = sum_sq_u + err_u^2;

            in_track = double(err_pos < pos_threshold);

            log_states(k,:) = [k, current_x', ref_now'];
            log_inputs(k,:) = [k, u_applied', uref_now'];

            log_errors(k,:) = [k, err_pos, err_att, err_vel, err_rate, err_thrust, e_yaw_raw, e_yaw_wrap, err_full, sum_sq_pos/k, sqrt(sum_sq_pos/k), sum_sq_full/k, sqrt(sum_sq_full/k), sum_sq_u/k, sqrt(sum_sq_u/k), in_track, calc_time, double(deadline_missed), tx_time, total_actuation_delay, state_age, double(qp_failed)];

            e_ps = abs(e);
            e_ps(angle_idx) = e_ps(angle_idx) * 180/pi;
            log_perst(k,:) = e_ps';

            current_x = next_x;
        end

        % Terminal sample at t = T_end_use. No additional MPC solve is made there.
        k_final = N_steps;

        [ref_final, uref_final] = refAt(x_ref_use, u_ref_use, k_final, x_hover, u_hover);

        e_final = current_x - ref_final;
        e_final_yaw_raw = e_final(9);
        e_final_yaw_wrap = atan2(sin(e_final_yaw_raw), cos(e_final_yaw_raw));
        e_final(9) = e_final_yaw_wrap;

        err_pos_final = norm(e_final(1:3));
        err_vel_final = norm(e_final(4:6));
        err_att_final = norm(e_final(7:9)) * 180/pi;
        err_rate_final = norm(e_final(10:12));
        err_thrust_final = norm(e_final(13:16));
        err_full_final = norm(e_final);

        sum_sq_pos_final = sum_sq_pos + err_pos_final^2;
        sum_sq_full_final = sum_sq_full + err_full_final^2;
        in_track_final = double(err_pos_final < pos_threshold);

        log_states(k_final,:) = [k_final, current_x', ref_final'];
        log_inputs(k_final,:) = [k_final, u_applied', uref_final'];

        log_errors(k_final,:) = [k_final, err_pos_final, err_att_final, err_vel_final, err_rate_final, err_thrust_final, e_final_yaw_raw, e_final_yaw_wrap, err_full_final, sum_sq_pos_final/N_steps, sqrt(sum_sq_pos_final/N_steps), sum_sq_full_final/N_steps, sqrt(sum_sq_full_final/N_steps), sum_sq_u/N_intervals, sqrt(sum_sq_u/N_intervals), in_track_final, NaN, NaN, NaN, NaN, NaN, NaN];

        e_ps = abs(e_final);
        e_ps(angle_idx) = e_ps(angle_idx) * 180/pi;
        log_perst(k_final,:) = e_ps';


        %% POST-LOOP ANALYSIS
        in_track_vec = logical(log_errors(1:N_steps_base, 16));

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
            rmse_ss_pos = sqrt(mean(log_errors(in_track_vec, 2).^2));
            rmse_ss_full = sqrt(mean(log_errors(in_track_vec, 9).^2));
        else
            rmse_ss_pos = NaN;
            rmse_ss_full = NaN;
        end

        max_pos_err = max(log_errors(1:N_steps_base, 2));
        max_yaw_err = max(abs(log_errors(1:N_steps_base, 8)));
        max_att_err = max(log_errors(1:N_steps_base, 3));

        solve_ref_times = solve_times(1:min(N_intervals_base, numel(solve_times)));
        solve_ref_times = solve_ref_times(~isnan(solve_ref_times));

        avg_solve = mean(solve_ref_times) * 1e3;
        worst_solve = max(solve_ref_times) * 1e3;
        std_solve = std(solve_ref_times) * 1e3;
        n_slow = sum(solve_ref_times > Ts * 0.8);

        n_deadline_miss = sum(deadline_miss(1:min(N_intervals_base, numel(deadline_miss))));

        deadline_miss_rate = n_deadline_miss / max(1, N_intervals_base);

        n_qp_fail = sum(qp_fail_log(1:min(N_intervals_base, numel(qp_fail_log))));

        tx_ref_times = tx_times(1:min(N_intervals_base, numel(tx_times)));
        total_lat_ref = total_latencies(1:min(N_intervals_base, numel(total_latencies)));

        avg_tx = mean(tx_ref_times, 'omitnan') * 1e3;
        avg_total_latency = mean(total_lat_ref, 'omitnan') * 1e3;

        perstate_rmse = zeros(1,16);

        for si = 1:16
            perstate_rmse(si) = sqrt(mean(log_perst(1:N_steps_base, si).^2));
        end

        rmse_pos_ref = log_errors(N_steps_base, 11);
        rmse_full_ref = log_errors(N_steps_base, 13);

        input_err_vec = vecnorm(log_inputs(1:N_intervals_base, 2:5) - log_inputs(1:N_intervals_base, 6:9), 2, 2);

        rmse_input_ref = sqrt(mean(input_err_vec.^2));
        mse_pos_ref = log_errors(N_steps_base, 10);
        mse_full_ref = log_errors(N_steps_base, 12);
        mse_input_ref = mean(input_err_vec.^2);

        fprintf(['  samples=%d (ref=%d)  intervals=%d  avg=%.1fms  worst=%.1fms  ', 'slow=%d  deadline_miss=%d (%.1f%%)  qp_fail=%d\n'], N_steps, N_steps_base, N_intervals, avg_solve, worst_solve, n_slow, n_deadline_miss, 100*deadline_miss_rate, n_qp_fail);

        fprintf('  tx=%.1fms  total_latency=%.1fms  state_age=%s\n', avg_tx, avg_total_latency, ternaryStateAge(perfect_state, state_acquisition_delay));

        fprintf('  RMSE_pos=%.4fm  RMSE_ss=%.4fm  t_conv=%.2fs \n\n', rmse_pos_ref, rmse_ss_pos, t_conv);

        if n_qp_fail > 0
            warning('C00_MPCSimulationLoop:QPFail', '%s | %s: %d unsolved QPs on the reference span (the run may have diverged).', scenario_name, cond_name, n_qp_fail);
        end


        %% FIGURES
        fig_title = sprintf('%s — %s', scenario_name, cond_name);

        %% Figure 1: State Tracking
        fig1 = figure('Name', sprintf('States: %s', fig_title), 'Color', 'w', 'Position', [50 50 1400 950]);

        for i = 1:16
            subplot(4,4,i);

            sc = 1;

            if ismember(i, angle_idx)
                sc = 180/pi;
            end

            plot(t_steps, log_states(:,i+1)*sc, 'b', 'LineWidth', 1.5);

            hold on;
            grid on;

            plot(t_steps, log_states(:,i+17)*sc, 'r--', 'LineWidth', 1.2);

            % End of tracking / beginning of hover
            xline(T_end_base, 'k:', 'LineWidth', 1.2);

            % Label only on the first subplot to avoid clutter
            if i == 1
                xline(T_end_base, 'k:', 'Tracking end / Hover start', 'LineWidth', 1.2, 'LabelVerticalAlignment', 'middle', 'LabelHorizontalAlignment', 'left', 'FontSize', 8);
            end

            ylabel(state_labels{i}, 'Interpreter', 'latex');

            xlim([0 T_end_use]);

            if i >= 13
                xlabel('t [s]');
            end

            if i == 1
                legend('Simulated', 'Reference', 'Location', 'best', 'Interpreter', 'none');
            end
        end

        sgtitle(sprintf('States Tracking: %s', fig_title), 'Interpreter', 'none');

        saveas(fig1, fullfile(run_dir, 'plot_states.png'));


        %% Figure 2: Applied Inputs -- FULL duration (tracking + hover tail)
        fig2 = figure('Name', sprintf('Inputs: %s', fig_title), 'Color', 'w', 'Position', [150 100 800 600]);

        for i = 1:4
            subplot(2,2,i);

            plot(t_steps, log_inputs(:,i+1), 'b', 'LineWidth', 1.5);

            hold on;
            grid on;

            plot(t_steps, log_inputs(:,i+5), 'r--', 'LineWidth', 1.2);

            % End of tracking / beginning of hover
            xline(T_end_base, 'k:', 'LineWidth', 1.2);

            % Label only on the first subplot
            if i == 1
                xline(T_end_base, 'k:', 'Tracking end / Hover start', 'LineWidth', 1.2, 'LabelVerticalAlignment', 'middle', 'LabelHorizontalAlignment', 'left', 'FontSize', 8);
            end

            % Input constraints
            yline(cp.u_min_16(i), 'k:', 'LineWidth', 0.8);

            yline(cp.u_max_16(i), 'k:', 'LineWidth', 0.8);

            % Hover thrust
            yline(T_hover_pm, 'g:', 'LineWidth', 0.8);

            ylim([cp.u_min_16(i), cp.u_max_16(i)]);
            ylabel(input_labels{i});

            xlim([0 T_end_use]);

            if i >= 3
                xlabel('t [s]');
            end

            if i == 1
                legend('Applied', 'Reference', 'Location', 'best', 'Interpreter', 'none');
            end
        end

        sgtitle(sprintf('Inputs Applied: %s', fig_title), 'Interpreter', 'none');

        saveas(fig2, fullfile(run_dir, 'plot_inputs.png'));


        %% Figure 3: Tracking Analysis -- REFERENCE SPAN ONLY
        %
        % Only the actual tracking interval [0, T_end_base] is shown.
        % The subsequent hover tail is excluded from the time-domain analysis.

        t_steps_ref = t_steps(1:N_steps_base);
        solve_ref = solve_times(1:N_intervals_base);

        fig3 = figure('Name', sprintf('Analysis: %s', fig_title), 'Color', 'w', 'Position', [250 100 1000 900]);


        %% 1. Position error
        subplot(3,2,1);

        semilogy(t_steps_ref, log_errors(1:N_steps_base,2) + 1e-12, 'b', 'LineWidth', 1.5);

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


        %% 2. Attitude error
        subplot(3,2,2);

        plot(t_steps_ref, log_errors(1:N_steps_base,3), 'Color', [0.8 0.2 0.2], 'LineWidth', 1.5);

        grid on;

        xlim([0 T_end_base]);

        ylabel('Attitude error [deg]');
        xlabel('t [s]');
        title('Attitude Error  (\phi, \theta, \psi)');


        %% 3. Yaw error
        subplot(3,2,3);

        plot(t_steps_ref, log_errors(1:N_steps_base,7)*180/pi, 'Color', [0.6 0.6 0.6], 'LineWidth', 1.0);

        hold on;
        grid on;

        plot(t_steps_ref, log_errors(1:N_steps_base,8)*180/pi, 'Color', [0.1 0.4 0.8], 'LineWidth', 1.8);

        xlim([0 T_end_base]);

        legend('algebraic', 'corrected', 'Location', 'best');

        ylabel('[deg]');
        xlabel('t [s]');
        title('Yaw Error');


        %% 4. Motor thrust error
        subplot(3,2,4);

        plot(t_steps_ref, log_errors(1:N_steps_base,6), 'Color', [0.5 0.2 0.6], 'LineWidth', 1.5);

        grid on;

        xlim([0 T_end_base]);

        ylabel('Thrust error [N]');
        xlabel('t [s]');
        title('Motor Thrust Error  (T_1..T_4 vs NLP reference)');


        %% 5. Per-state RMSE
        subplot(3,2,5);

        hold on;
        grid on;

        h1 = bar(1:3, perstate_rmse(1:3), 'FaceColor', [0.25 0.45 0.80]);

        h2 = bar(4:6, perstate_rmse(4:6), 'FaceColor', [0.20 0.65 0.30]);

        h3 = bar(7:9, perstate_rmse(7:9), 'FaceColor', [0.85 0.50 0.10]);

        h4 = bar(10:12, perstate_rmse(10:12), 'FaceColor', [0.75 0.20 0.20]);

        h5 = bar(13:16, perstate_rmse(13:16), 'FaceColor', [0.55 0.25 0.65]);

        xticks(1:16);

        ax5 = gca;
        ax5.TickLabelInterpreter = 'latex';

        ax5.XTickLabel = {'$r_n$', '$r_e$', '$r_d$', '$v_n$', '$v_e$', '$v_d$', '$\phi$', '$\theta$', '$\psi$', '$p$', '$q$', '$r$', '$T_1$','$T_2$','$T_3$','$T_4$'};

        xtickangle(30);

        ylabel('RMSE');
        xlabel('State');

        title('Per-State RMSE (tracking span)');

        legend([h1 h2 h3 h4 h5], {'position','velocity','angles','rates','thrust'}, 'Location', 'northeast', 'FontSize', 7);


        %% 6. Controller solve time
        subplot(3,2,6);

        plot(t_steps_ref(1:end-1), solve_ref*1e3, 'Color', [0.3 0.3 0.3], 'LineWidth', 0.8);

        hold on;
        grid on;

        yline(Ts*1e3, 'r--', sprintf('T_s = %.0f ms', Ts*1e3), 'LabelHorizontalAlignment', 'left');

        yline(Ts*0.8*1e3, 'r:', '80%% budget', 'LabelHorizontalAlignment', 'right');

        xlim([0 T_end_base]);

        ylabel('Solve time [ms]');
        xlabel('t [s]');

        title(sprintf('Controller time  (worst = %.1f ms, tracking span)', max(solve_ref)*1e3));


        %% Overall title
        sgtitle(sprintf('Tracking Analysis: %s (tracking span only, hover excluded)', fig_title), 'Interpreter', 'none');

        saveas(fig3, fullfile(run_dir, 'plot_analysis.png'));


        %% GIF
        if make_gif
            gif_name = fullfile(run_dir, 'tracking.gif');

            GenerateTrackingSimulationGIF(log_states, x_ref_use, gif_name, fig_title, '3D', mp, C00_PlotDrone);
        end


        %% CSV, full duration logs kept complete
        writetable(array2table(log_states, 'VariableNames', state_headers), fullfile(run_dir, 'states.csv'));

        writetable(array2table(log_inputs, 'VariableNames', input_headers), fullfile(run_dir, 'inputs.csv'));

        writetable(array2table(log_errors, 'VariableNames', error_headers), fullfile(run_dir, 'errors.csv'));

        ps_names = arrayfun(@(i) sprintf('state_%d_error', i), 1:16, 'UniformOutput', false);

        writetable(array2table(log_perst, 'VariableNames', ps_names), fullfile(run_dir, 'per_state_errors.csv'));


        %% ACCUMULATE METRICS, reference span only
        metrics_summary(end+1,:) = {cond_name, scenario_name, rmse_pos_ref, rmse_full_ref, rmse_input_ref, avg_solve, N_steps_base, T_end_base};  %#ok<SAGROW>

        metrics_extended(end+1,:) = {cond_name, scenario_name, rmse_pos_ref, rmse_ss_pos, rmse_full_ref, rmse_ss_full, mse_pos_ref, mse_full_ref, t_conv, max_pos_err, max_yaw_err*180/pi, max_att_err, rmse_input_ref, mse_input_ref, avg_solve, worst_solve, std_solve, n_slow, n_deadline_miss, deadline_miss_rate, avg_tx, avg_total_latency, n_qp_fail, N_steps_base, T_end_base};  %#ok<SAGROW>

    end
end


%% EXPORT SUMMARY
cols_sum = {'initial_condition','scenario', 'RMSE_position_m','RMSE_fullstate','RMSE_input', 'avg_solve_time_ms','steps_reference_span','reference_duration_s'};

writetable(cell2table(metrics_summary, 'VariableNames', cols_sum), fullfile(fullRunDir, 'MPC_metrics_summary.csv'));


%% EXPORT EXTENDED
cols_ext = {'initial_condition','scenario', 'rmse_position_all_steps_m','rmse_position_steadystate_m', 'rmse_fullstate_all_steps','rmse_fullstate_steadystate', 'mse_position_m2','mse_fullstate', 'convergence_time_s', 'peak_position_error_m','peak_yaw_error_deg','peak_attitude_error_deg', 'rmse_input','mse_input', 'avg_solve_time_ms','worst_solve_time_ms','std_solve_time_ms', 'steps_exceeding_80pct_budget','deadline_misses','deadline_miss_rate', 'avg_tx_time_ms','avg_total_latency_ms','qp_failures', 'steps_reference_span','reference_duration_s'};

writetable(cell2table(metrics_extended, 'VariableNames', cols_ext), fullfile(fullRunDir, 'MPC_metrics_extended.csv'));

fprintf('Done. Results: %s\n', fullRunDir);


%% DEBUG: Omega_i plausibility check across all (IC x scenario) combinations
% Converts logged T1..T4 [N] to Omega1..Omega4 [rad/s] via Omega = sqrt(T/K_T).
K_T = mp.K_T;
Omega_min = sqrt(cp.u_min_16(1) / K_T);
Omega_max = sqrt(cp.u_max_16(1) / K_T);

motor_colors = {[0.20 0.45 0.80], [0.85 0.20 0.20], [0.20 0.65 0.30], [0.60 0.30 0.70]};

figure('Name', 'Omega_i plausibility check (all combos)', 'Color', 'w', 'Position', [50 50 1500 1100]);
tl = tiledlayout(numel(initial_conditions), size(scenarios, 1), 'Padding', 'compact', 'TileSpacing', 'compact');

for c = 1:numel(initial_conditions)
    cond_name = initial_conditions{c};
    for s = 1:size(scenarios, 1)
        scenario_name = scenarios{s, 1};

        run_dir = fullfile(fullRunDir, scenario_name, cond_name);
        Tstates = readtable(fullfile(run_dir, 'states.csv'));
        t_ = (Tstates.step - 1) * Ts;

        nexttile;
        plot(t_, sqrt(max(Tstates.T1, 0) / K_T), 'Color', motor_colors{1}, 'LineWidth', 1.1); hold on; grid on;
        plot(t_, sqrt(max(Tstates.T2, 0) / K_T), 'Color', motor_colors{2}, 'LineWidth', 1.1);
        plot(t_, sqrt(max(Tstates.T3, 0) / K_T), 'Color', motor_colors{3}, 'LineWidth', 1.1);
        plot(t_, sqrt(max(Tstates.T4, 0) / K_T), 'Color', motor_colors{4}, 'LineWidth', 1.1);

        yline(Omega_min, 'k:', 'LineWidth', 0.7);
        yline(Omega_max, 'k:', '$\Omega_{\max}$', 'Interpreter', 'latex', 'LineWidth', 0.7, 'LabelHorizontalAlignment', 'right');

        yline(sqrt(T_hover_pm / K_T), 'g:', '$\Omega_{\mathrm{hover}}$', 'Interpreter', 'latex', 'LineWidth', 0.7, 'LabelHorizontalAlignment', 'right');

        title(sprintf('%s — %s', scenario_name, cond_name), 'FontSize', 9);

        if s == 1
            ylabel('$\Omega$ [rad/s]', 'Interpreter', 'latex');
        end

        if c == numel(initial_conditions)
            xlabel('t [s]');
        end
    end
end

lgd = legend({'$\Omega_1$', '$\Omega_2$', '$\Omega_3$', '$\Omega_4$'}, 'Interpreter', 'latex', 'Orientation', 'horizontal');

lgd.Layout.Tile = 'south';

sgtitle(tl, '$\mathrm{Motor\ angular\ velocities:}\quad \Omega_i = \sqrt{T_i/K_T}$', 'Interpreter', 'latex');


%% LOCAL FUNCTIONS
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


function x_meas = delayedState(t_hist, x_hist, t_query)
    % DELAYEDSTATE  Linear interpolation of the plant state at an older time
    %               stamp, on the DENSE ode45 history (t_hist: 1 x n, x_hist: 16 x n).
    %               The query is clamped to the available simulation history.
    t_query = min(max(t_query, t_hist(1)), t_hist(end));

    if numel(t_hist) == 1
        x_meas = x_hist(:,1);
        return;
    end

    x_meas = interp1(t_hist(:), x_hist', t_query, 'linear')';
end


function out = ternaryStateAge(perfect_state, delayed_value)
    % TERNARYSTATEAGE  Compact text helper for console reporting.
    if perfect_state
        out = '0 ms';
    else
        out = sprintf('%.1f ms', delayed_value*1e3);
    end
end