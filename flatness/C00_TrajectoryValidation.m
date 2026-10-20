function report = C00_TrajectoryValidation(x_ref, u_ref, t_steps, info, name, do_plot)
% Validates an NLP-generated reference before it is used by the MPC.
%
% Question answered: if the four motors are driven by u_ref (PWC, ZOH),
%   1. do they stay inside their physical limits, with room left for the MPC?
%   2. does the drone really follow the reference, also BETWEEN the nodes?
%   3. does every state stay inside the MPC state box?
%
% The reference is iterated over the WHOLE trajectory: the PWC input is
% replayed open loop with ode45, evaluated on the fine grid of the continuous reference
% (M points per interval), and compared with x^ref(t).
%
% ACCEPTANCE THRESHOLDS = the tolerances (sigma) used as Bryson weights in
% the NLP (info.opts.sig_x, sig_du, sig_ddu). The weights say "an error of
% one sigma costs 1"; this script says "an error larger than sigma FAILS".
%
% INPUTS
%   x_ref   - nx x (N+1)  optimized states at the nodes
%   u_ref   -  4 x (N+1)  optimized input (last column is padding, not applied)
%   t_steps -  1 x (N+1)  node times [s]
%   info    - struct returned by C00_GeneratePWC_reference
%             (uses xr_dense, opts, du, stats, max_defect)
%   name    - (optional) label for messages/figures
%   do_plot - (optional) logical, default true
%
% OUTPUT
%   report  - struct:
%       pass          true if no check is FAIL (WARN is allowed)
%       n_fail, n_warn
%       checks        cell array {name, value, limit, status}
%       ratio         struct with max |error|/sigma per state group
%       headroom      struct: up_min, dn_min [N] to the physical limits
%       authority     struct: up_min, dn_min (4x1) worst-case room left for
%                     the MPC on [collective m/s^2; roll; pitch; yaw rad/s^2]
%       x_rep, t_dense  ode45 replay on the fine grid

    if nargin < 5 || isempty(name),    name = 'trajectory'; end
    if nargin < 6 || isempty(do_plot), do_plot = true;      end

    %% Setup
    mp    = ModelParameters();
    model = C00_ANTX_quadcopter(mp, info.opts.delay_motors);

    nx = size(x_ref, 1);
    N  = numel(t_steps) - 1;
    Ts = t_steps(2) - t_steps(1);
    M  = info.opts.M;
    h  = Ts / M;

    U       = u_ref(:, 1:N);              % applied inputs (drop padded last column)
    sig_x   = info.opts.sig_x(1:nx);      % nx x 1, same tolerances as the NLP weights
    sig_du  = info.opts.sig_du;           % [N]
    sig_ddu = info.opts.sig_ddu;          % [N/s]

    u_min = ControlParameters.u_min_16;   u_max = ControlParameters.u_max_16;
    span  = u_max - u_min;
    u_lo  = u_min + info.opts.margin_frac * span;     % tightened bounds used by the NLP
    u_hi  = u_max - info.opts.margin_frac * span;
    x_min = ControlParameters.x_min_16(1:nx);
    x_max = ControlParameters.x_max_16(1:nx);
    T_hover = mp.m * mp.g / 4;

    state_names = {'r_n','r_e','r_d','v_n','v_e','v_d','phi','theta','psi', ...
                   'p','q','r','T1','T2','T3','T4'};
    state_names = state_names(1:nx);

    groups = {'position', 1:3; 'velocity', 4:6; 'attitude', 7:9; ...
              'body rates', 10:12; 'motor thrust states', 13:16};
    groups = groups(cellfun(@(ix) max(ix) <= nx, groups(:, 2)), :);


    %% Open-loop replay of u_ref on the fine grid (ode45)
    %
    % Interval k covers fine-grid columns (k-1)*M+1 ... k*M+1. Same grid as
    % info.xr_dense / info.t_dense.

    ode_opts = odeset('RelTol', 1e-9, 'AbsTol', 1e-10);
    x_rep = zeros(nx, N*M + 1);
    x_rep(:, 1) = x_ref(:, 1);
    x = x_ref(:, 1);

    for k = 1:N
        tk = t_steps(k) + (0:M) * h;                       % M+1 instants
        [~, xs] = ode45(@(t, x) model.dynamics(t, x, U(:, k)), tk, x, ode_opts);
        x_rep(:, (k-1)*M+1 : k*M+1) = xs.';
        x = xs(end, :).';
    end


    %% Errors, normalized by the tolerances (sigma)

    E  = x_rep - info.xr_dense(1:nx, :);      % replay vs CONTINUOUS reference (fine grid)
    Rn = abs(E) ./ sig_x;                     % 1.0 = error equal to its tolerance

    ratio = zeros(size(groups, 1), 1);
    for g = 1:size(groups, 1)
        ratio(g) = max(max(Rn(groups{g, 2}, :)));
    end

    En         = x_rep(:, 1:M:end) - x_ref;   % replay vs NLP nodes (consistency of the NLP)
    ratio_cons = max(max(abs(En) ./ sig_x));


    %% Thrust: headroom to the physical limits and MPC authority

    hr_up = u_max - U;                        % room to increase each motor  [N]
    hr_dn = U - u_min;                        % room to decrease each motor  [N]

    % Virtual inputs [T; L; M; N] = A * [T1..T4] 
    a = mp.b / sqrt(2);
    A = [ 1,  1,  1,  1;
          a*[ 1, -1, -1,  1];
          a*[ 1,  1, -1, -1];
          mp.Beta*[-1,  1, -1,  1] ];
    inertia = [mp.m; mp.Jxx; mp.Jyy; mp.Jzz];         % -> [a_z; p_dot; q_dot; r_dot]

    D_hi = u_max - U;   D_lo = u_min - U;
    auth_up = zeros(4, N);   auth_dn = zeros(4, N);
    for i = 1:4
        Pu = A(i, :).' .* D_hi;
        Pd = A(i, :).' .* D_lo;
        auth_up(i, :) = sum(max(Pu, Pd), 1) / inertia(i);     % max increase of channel i alone
        auth_dn(i, :) = sum(min(Pu, Pd), 1) / inertia(i);     % max decrease of channel i alone
    end
    auth_up_min = min(auth_up, [], 2);
    auth_dn_min = min(-auth_dn, [], 2);


    %% Checks

    slack = min(x_rep - x_min, x_max - x_rep);         % >= 0 inside the MPC box
    [worst_slack, i_w] = min(min(slack, [], 2));

    du = info.du;                                      % 4 x N correction
    r_du = max(abs(du(:))) / sig_du;

    checks = {};
    checks(end+1, :) = {'IPOPT converged', info.stats.return_status, 'Solve_Succeeded', ...
                        status(logical(info.stats.success), 'FAIL')};
    checks(end+1, :) = {'Shooting defect (max)', sprintf('%.1e', info.max_defect), '<= 1e-6', ...
                        status(info.max_defect <= 1e-6, 'FAIL')};
    checks(end+1, :) = {'Motor thrust inside PHYSICAL limits', ...
                        sprintf('[%.3f, %.3f] N', min(U(:)), max(U(:))), ...
                        sprintf('[%.3f, %.3f] N', u_min(1), u_max(1)), ...
                        status(all(all(U >= u_min & U <= u_max)), 'FAIL')};
    checks(end+1, :) = {'Motor thrust inside TIGHTENED bounds', ...
                        sprintf('[%.3f, %.3f] N', min(U(:)), max(U(:))), ...
                        sprintf('[%.3f, %.3f] N', u_lo(1), u_hi(1)), ...
                        status(all(all(U >= u_lo & U <= u_hi)), 'WARN')};
    checks(end+1, :) = {'States inside MPC box (fine grid)', ...
                        sprintf('tightest: %s (slack %.3g)', state_names{i_w}, worst_slack), ...
                        'slack >= 0', status(worst_slack >= -1e-9, 'FAIL')};
    checks(end+1, :) = {'NLP nodes vs ode45 replay', sprintf('%.3f sigma', ratio_cons), ...
                        '<= 0.10 sigma', status(ratio_cons <= 0.10, 'FAIL')};
    for g = 1:size(groups, 1)
        checks(end+1, :) = {['Tracking of x^r(t): ' groups{g, 1}], ...
                            sprintf('%.3f sigma', ratio(g)), '<= 1 sigma', ...
                            status(ratio(g) <= 1, 'FAIL')};
    end
    checks(end+1, :) = {'Correction size max|du|', ...
                        sprintf('%.2f sigma_du (%.2f mN)', r_du, 1e3*max(abs(du(:)))), ...
                        '<= 1 sigma_du', status(r_du <= 1, 'WARN')};
    if N > 1
        r_ddu = max(max(abs(diff(du, 1, 2)))) / (sig_ddu * Ts);
        checks(end+1, :) = {'Correction smoothness max|du_k+1 - du_k|', ...
                            sprintf('%.2f sigma_ddu*Ts', r_ddu), '<= 1', ...
                            status(r_ddu <= 1, 'WARN')};
    end

    n_fail = sum(strcmp(checks(:, 4), 'FAIL'));
    n_warn = sum(strcmp(checks(:, 4), 'WARN'));


    %% Report

    fprintf('\n Trajectory validation: %s (N = %d, Ts = %.3f s)\n', name, N, Ts);
    for i = 1:size(checks, 1)
        fprintf('  [%s] %-44s %-34s (limit: %s)\n', checks{i, 4}, checks{i, 1}, checks{i, 2}, checks{i, 3});
    end

    fprintf('  Headroom to physical limits: up %.3f N (%.0f%% of hover) | down %.3f N (%.0f%% of hover)\n', ...
        min(hr_up(:)), 100*min(hr_up(:))/T_hover, min(hr_dn(:)), 100*min(hr_dn(:))/T_hover);
    fprintf('  Peak collective thrust: %.3f N = %.2f x weight (available: %.2f x weight)\n', ...
        max(sum(U, 1)), max(sum(U, 1))/(mp.m*mp.g), sum(u_max)/(mp.m*mp.g));
    fprintf('  Worst-case authority left for the MPC (each channel alone), up / down:\n');
    ch_names = {'collective [m/s^2]', 'roll  acc [rad/s^2]', 'pitch acc [rad/s^2]', 'yaw   acc [rad/s^2]'};
    for i = 1:4
        fprintf('    %-22s +%9.2f / -%9.2f\n', ch_names{i}, auth_up_min(i), auth_dn_min(i));
    end
    if n_fail == 0
        fprintf('  VERDICT: PASS (%d warning(s))\n\n', n_warn);
    else
        fprintf('  VERDICT: FAIL (%d failed check(s), %d warning(s))\n\n', n_fail, n_warn);
    end


    %% Output struct

    report.pass      = (n_fail == 0);
    report.n_fail    = n_fail;
    report.n_warn    = n_warn;
    report.checks    = checks;
    report.ratio     = cell2struct(num2cell(ratio), regexprep(groups(:, 1), '\s', '_'), 1);
    report.headroom  = struct('up_min', min(hr_up(:)), 'dn_min', min(hr_dn(:)));
    report.authority = struct('up_min', auth_up_min, 'dn_min', auth_dn_min);
    report.x_rep     = x_rep;
    report.t_dense   = info.t_dense;


    %% Figure

    if do_plot
        figure('Name', sprintf('%s Validation (C00)', name), 'Color', 'w', 'Position', [30 50 1000 650]);

        subplot(2, 1, 1);
        cols = lines(4);
        for j = 1:4
            stairs(t_steps(1:N), U(j, :), 'Color', cols(j, :), 'LineWidth', 1.4); hold on;
        end
        yline(u_min(1), 'r-');  yline(u_max(1), 'r-',  'physical');
        yline(u_lo(1),  'r--'); yline(u_hi(1),  'r--', 'tightened');
        yline(T_hover,  'k:',  'T_{hover}/4');
        ylim([u_min(1), 1.05*u_max(1)]);
        ylabel('u_{ref} [N]'); grid on;
        legend('T_1', 'T_2', 'T_3', 'T_4', 'Location', 'east');
        title('Motor thrust vs limits');

        subplot(2, 1, 2);
        for g = 1:size(groups, 1)
            semilogy(info.t_dense, max(Rn(groups{g, 2}, :), [], 1), 'LineWidth', 1.3); hold on;
        end
        yline(1, 'r--', '\sigma (tolerance)');
        ylabel('max |error| / \sigma'); xlabel('t [s]'); grid on;
        legend(groups(:, 1)', 'Location', 'best');
        title('Tracking of the continuous reference (ode45 replay, fine grid)');

        sgtitle(sprintf('%s: reference validation | %s', name, verdict_str(report.pass, n_warn)));
    end

end


function st = status(is_ok, level)
% STATUS  'PASS' if the condition holds, otherwise the given level ('FAIL'|'WARN').
    if is_ok
        st = 'PASS';
    else
        st = level;
    end
end


function s = verdict_str(pass, n_warn)
% VERDICT_STR  Short text for figure titles.
    if pass
        s = sprintf('PASS (%d warning(s))', n_warn);
    else
        s = 'FAIL';
    end
end
