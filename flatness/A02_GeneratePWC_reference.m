function [x_ref, u_ref, t_steps, info] = A02_GeneratePWC_reference(traj_obj, T_sim, Ts, nlp_opts, quad_options)
% Generates an OPTIMIZED Piece-Wise Constant (PWC) reference
% for the A02 model, by solving a nonlinear program (NLP) that
% corrects the interval-averaged input: u_k = u_bar_k + du_k*
%
%     Flatness map : A02_FlatnessMap (continuous reference x^r(t))
%     Model        : A02_SimplifiedModel_Deployment(mp)
%                    (12 states, inertial NED velocities, u = [T; L; M; N],
%                    same model used by the MPC)
%     Solver       : CasADi + IPOPT, multiple shooting, fixed-step RK4
%
% u_bar_k is the TIME-AVERAGE of the flatness input over [t_k, t_k+Ts]
% (computed via the auxiliary ODE dz/dt = u(t)).
% It is used as the starting point; the NLP finds the correction du_k.
%
% NLP (multiple shooting):
%
%   decision variables  w = { x_0 ... x_N , du_0 ... du_{N-1} }
%
%   min   J = sum_k sum_j h * (x_kj - x^r(t_kj))' Q (x_kj - x^r(t_kj))   tracking, also INSIDE intervals
%           + Ts * sum_k du_k' R_du du_k                                  keep correction small
%           + (1/Ts) * sum_k (du_{k+1}-du_k)' S (du_{k+1}-du_k)           no chattering
%
%   s.t.  x_0     = x^r(0)                                  initial condition
%         x_{k+1} = Phi(x_k, u_bar_k + du_k),  k=0..N-1      continuity (RK4, M substeps, ZOH)
%         x_N     = x^r(t_0)                                terminal / periodicity (if periodic)
%         T_lo   <= A^{-1}*(u_bar_k + du_k) <= T_hi            motor polytope with margin
%         x_min  <= x_k <= x_max                             state box (nodes only)
%
% x_kj are the RK4 substep states of interval k (j = 0..M-1).
% A^{-1} = inverse allocation [T;L;M;N] -> [T1..T4] (configuration C): the
% physical limits are on the single rotors, not on T, L, M, N.
%
% TUNING RULE (all weights are Bryson weights 1/sigma^2, sigma = tolerated
% deviation). ONE design number, ONE time scale, the rest is physics:
%   sig_p   : position tolerance [m]                      (design number - 0.01m)
%   tau_c   : error accumulation time [s], default 2*Ts   (time scale)
%   sig_v   = sig_p / tau_c                               velocity [m/s]
%   sig_th  = sig_v / (g*tau_c)                           tilt     [rad]
%   sig_w   = sig_th / tau_c                              rates    [rad/s]
%   sig_u   = [ m*sig_v/tau_c;                            T [N]     (collective -> velocity)
%               Jxx*sig_w/tau_c;                          L [N*m]   (roll  -> rate)
%               Jyy*sig_w/tau_c;                          M [N*m]   (pitch -> rate)
%               Jzz*sig_w/tau_c ]                         N [N*m]   (yaw   -> rate)
%   sig_du  = kappa_du * sig_u                            correction, per channel
%   sig_ddu = sig_du / tau_c                              rate of change of the correction
%
%
% INPUTS
%   traj_obj     - trajectory object (ShapeLemniscate | ShapeCircle | ...)
%                  must provide get_flat_outputs(t)
%   T_sim        - total duration [s]; should be an integer multiple of Ts
%   Ts           - sampling time [s] (0.1 for GCS deployment)
%   nlp_opts     - (optional) struct overriding any field of
%                  defaultOptions() below (weights, M, margin, ...)
%   quad_options - (optional) odeset for the u_bar quadrature ODE.
%                  Default: odeset('RelTol',1e-10,'AbsTol',1e-12)
%
% OUTPUTS
%   x_ref   - 12 x (N+1) optimized state trajectory at the nodes
%   u_ref   -  4 x (N+1) optimized input u_bar + du* [T;L;M;N]
%             (last column repeated, as in previous versions)
%   t_steps -  1 x (N+1) time vector [s], exact multiples of Ts
%   info    - struct with diagnostics:
%               u_bar      4 x (N+1) averaged input (v2, case B)
%               du         4 x N     optimal correction
%               J          struct: total, track, du, ddu
%               stats      IPOPT stats (success, return_status, iter_count)
%               max_defect max continuity violation (should be ~0)
%               t_dense, xr_dense  continuous reference at RK4 substeps
%               Ainv       4 x 4 inverse allocation
%               opts       options actually used (sig_x, sig_u, sig_du, sig_ddu, ...)
%
% REQUIREMENTS
%   CasADi (MATLAB build, includes IPOPT) on the MATLAB path.
%   A02_SimplifiedModel_Deployment.dynamics must accept CasADi symbols

    import casadi.*

    %% Optional arguments
    if nargin < 4 || isempty(nlp_opts)
        nlp_opts = struct();
    end
    if nargin < 5 || isempty(quad_options)
        quad_options = odeset('RelTol', 1e-10, 'AbsTol', 1e-12);
    end
    opts = mergeOptions(defaultOptions(), nlp_opts);

    %% Model parameters and model instance
    mp    = ModelParameters();
    model = A02_SimplifiedModel_Deployment(mp);

    nx = 12;
    nu = 4;

    %% Inverse allocation [T;L;M;N] -> [T1;T2;T3;T4] (configuration C)
    Ainv = allocationInverse(mp);

    %% Tolerances (Bryson sigmas) derived from ONE design tolerance.
    if isempty(opts.tau_c),   opts.tau_c   = 2*Ts;  end
    sg = deriveSigmas(mp, opts.sig_p, opts.tau_c, opts.kappa_du);
    if isempty(opts.sig_x),   opts.sig_x   = sg.x;   end
    if isempty(opts.sig_u),   opts.sig_u   = sg.u;   end
    if isempty(opts.sig_du),  opts.sig_du  = sg.du;  end
    if isempty(opts.sig_ddu), opts.sig_ddu = sg.ddu; end

    %% Time grid (exact multiples of Ts, as required by the GCS)
    N       = round(T_sim / Ts);     % number of shooting intervals
    t_steps = (0:N) * Ts;            % 1 x (N+1) nodes
    M       = opts.M;                % RK4 substeps per interval
    h       = Ts / M;                % RK4 step [s]

    if abs(N*Ts - T_sim) > 1e-9
        warning('A02_GeneratePWC_reference:Grid', ...
            'T_sim is not a multiple of Ts: grid ends at %.4f s instead of %.4f s.', N*Ts, T_sim);
    end


    %% Continuous reference x^r(t) on the dense RK4 grid
    %
    % t_sub is ordered interval by interval:
    %   [t_0, t_0+h, ..., t_0+(M-1)h, t_1, t_1+h, ...]
    % so that block k (M columns) holds the reference of interval k.

    t_sub   = reshape(t_steps(1:N) + (0:M-1)' * h, 1, []);   % 1 x N*M
    t_dense = [t_sub, t_steps(N+1)];                         % 1 x (N*M+1)

    xr_dense = zeros(nx, numel(t_dense));
    for i = 1:numel(t_dense)
        [s, ds, dds, ddds, dddds] = traj_obj.get_flat_outputs(t_dense(i));
        xr_dense(:, i) = A02_FlatnessMap.map(mp, s, ds, dds, ddds, dddds);
    end
    xr_dense(9, :) = unwrap(xr_dense(9, :));   % continuous yaw (no 2*pi jumps in the cost)

    XR      = xr_dense(:, 1:N*M);                     % nx x N*M  (interval blocks)
    xr_node = xr_dense(:, [1:M:N*M, N*M+1]);          % nx x (N+1) reference at nodes
    x0      = xr_node(:, 1);                          % initial condition
    xT      = xr_node(:, end);                        % terminal target


    %% Interval-averaged input u_bar (v2 reference)
    u_bar = zeros(nu, N);
    for k = 1:N
        u_bar(:, k) = computeMeanInput_ode45(traj_obj, mp, t_steps(k), Ts, quad_options);
    end


    %% Cost weights (Bryson: 1 / max acceptable deviation^2)

    q_diag = 1 ./ opts.sig_x(:).^2;     % nx x 1  (Q = diag(q_diag))
    r_du   = 1 ./ opts.sig_du(:).^2;    % nu x 1  (R_du = diag(r_du)), per channel
    s_ddu  = 1 ./ opts.sig_ddu(:).^2;   % nu x 1  (S = diag(s_ddu)),   per channel


    %% One-interval integrator + running cost (symbolic, SX)
    %
    % F_int(x_k, u_k, XR_k) -> [x_{k+1}, J_k]
    % mapped over the N intervals (evaluated all at once).

    F_int = buildIntervalFunction(model, nx, nu, M, h, q_diag);
    F_map = F_int.map(N);


    %% NLP (multiple shooting)

    X  = MX.sym('X',  nx, N+1);    % states at the nodes (decision variables)
    DU = MX.sym('DU', nu, N);      % input corrections   (decision variables)

    U = u_bar + DU;                              % applied ZOH input
    [X_next, J_int] = F_map(X(:, 1:N), U, XR);   % propagation of every interval

    % Cost
    J_track = sum2(J_int);
    J_du    = Ts * sum2(sum1(repmat(r_du, 1, N) .* DU.^2));
    dDU     = DU(:, 2:N) - DU(:, 1:N-1);
    J_ddu   = (1/Ts) * sum2(sum1(repmat(s_ddu, 1, N-1) .* dDU.^2));   % sig_ddu per channel
    J       = J_track + J_du + J_ddu;

    % Continuity constraints (shooting gaps = 0) + motor polytope (linear in U)
    g_cont = vec(X(:, 2:N+1) - X_next);
    g_poly = vec(Ainv * U);                      % [T1..T4] of every interval
    g      = [g_cont; g_poly];

    w   = [vec(X); vec(DU)];
    nlp = struct('x', w, 'f', J, 'g', g);


    %% Bounds

    % States: same box as the MPC; x_0 fixed; x_N fixed if periodic
    x_min = ControlParameters.x_min_12;
    x_max = ControlParameters.x_max_12;
    lbX = repmat(x_min, 1, N+1);
    ubX = repmat(x_max, 1, N+1);
    lbX(:, 1) = x0;   ubX(:, 1) = x0;
    if opts.periodic
        lbX(:, end) = xT;   ubX(:, end) = xT;
    end

    % Corrections: free (the physical limit is the motor polytope below)
    lbDU = -inf(nu, N);
    ubDU =  inf(nu, N);

    % Motor polytope: single-rotor limits tightened by a margin (leave authority to the MPC)
    span = ControlParameters.u_max_16 - ControlParameters.u_min_16;
    T_lo = ControlParameters.u_min_16 + opts.margin_frac * span;     % 4 x 1 [N]
    T_hi = ControlParameters.u_max_16 - opts.margin_frac * span;

    Tm_bar = Ainv * u_bar;                       % rotor thrusts of u_bar
    if any(any(Tm_bar < T_lo)) || any(any(Tm_bar > T_hi))
        warning('A02_GeneratePWC_reference:BarOutsideBounds', ...
            'u_bar violates the tightened motor bounds at some steps (du = 0 not feasible there).');
    end

    lbw = [lbX(:); lbDU(:)];
    ubw = [ubX(:); ubDU(:)];
    lbg = [zeros(nx*N, 1); repmat(T_lo, N, 1)];
    ubg = [zeros(nx*N, 1); repmat(T_hi, N, 1)];


    %% Initial guess and solve
    %
    % States on the IDEAL continuous reference, corrections at zero (i.e. start from u_bar).

    w0 = [xr_node(:); zeros(nu*N, 1)];

    s_opts = struct();
    s_opts.ipopt.max_iter    = opts.max_iter;
    s_opts.ipopt.tol         = opts.tol;
    s_opts.ipopt.print_level = opts.print_level;
    s_opts.print_time        = opts.print_level > 0;

    solver = nlpsol('solver', 'ipopt', nlp, s_opts);
    sol    = solver('x0', w0, 'lbx', lbw, 'ubx', ubw, 'lbg', lbg, 'ubg', ubg);
    stats  = solver.stats();

    if ~stats.success
        warning('A02_GeneratePWC_reference:NotConverged', ...
            'IPOPT did not converge: %s', stats.return_status);
    end


    %% Extract solution

    w_opt  = full(sol.x);
    X_opt  = reshape(w_opt(1:nx*(N+1)), nx, N+1);
    DU_opt = reshape(w_opt(nx*(N+1)+1:end), nu, N);

    x_ref = X_opt;
    u_ref = [u_bar + DU_opt, u_bar(:, end) + DU_opt(:, end)];   % repeat last input

    % Cost breakdown (diagnostics)
    J_parts = Function('J_parts', {w}, {J_track, J_du, J_ddu});
    [jt, jd, jdd] = J_parts(w_opt);

    g_opt = full(sol.g);

    info.u_bar      = [u_bar, u_bar(:, end)];
    info.du         = DU_opt;
    info.J          = struct('total', full(sol.f), 'track', full(jt), ...
                             'du', full(jd), 'ddu', full(jdd));
    info.stats      = stats;
    info.max_defect = max(abs(g_opt(1:nx*N)));   % continuity part only
    info.t_dense    = t_dense;
    info.xr_dense   = xr_dense;
    info.Ainv       = Ainv;
    info.opts       = opts;

    fprintf('\n[A02_GeneratePWC_reference] %s | iter %d | J = %.4g (track %.4g, du %.4g, ddu %.4g)\n', ...
        stats.return_status, stats.iter_count, info.J.total, info.J.track, info.J.du, info.J.ddu);
    fprintf('  max|du| = [%.4f N, %.2e %.2e %.2e N*m] | max defect = %.2e | N = %d intervals, M = %d substeps\n', ...
        max(abs(DU_opt), [], 2), info.max_defect, N, M);
    fprintf('  sigma: pos %.3g m | u [%.3g N, %.2e %.2e %.2e N*m] | du = %g*sigma_u  (tau_c = %.2f s)\n\n', ...
        opts.sig_x(1), opts.sig_u, opts.kappa_du, opts.tau_c);

end


function opts = defaultOptions()
% DEFAULTOPTIONS  Default NLP settings. Any field can be overridden via
%                 the nlp_opts input.
%
% Weights follow Bryson's rule: weight = 1 / sigma^2, where sigma is the
% largest deviation considered acceptable for that quantity.

    opts.M            = 10;      % RK4 substeps per interval (h = Ts/M)
    opts.periodic     = true;    % x_N = x^r(t_N); set false for open curves (spiral)
    opts.margin_frac  = 0.10;    % motor margin, fraction of [u_min, u_max] span (per rotor)

    % Tolerances (Bryson sigmas). See TUNING RULE in the header.
    opts.sig_p    = 0.01;   % [m]  the ONE design tolerance (position error of the reference)
    opts.tau_c    = [];     % [s]  error accumulation time scale; [] -> 2*Ts
    opts.kappa_du = 10;     % [-]  sig_du = kappa_du * sig_u  (correction is cheap: regularizer only)

    % Overrides ([] = derived by the rule)
    opts.sig_x   = [];      % 12 x 1, state order of A02_SimplifiedModel_Deployment
    opts.sig_u   = [];      % 4 x 1  [N; N*m; N*m; N*m] input tolerance (used by the MPC R)
    opts.sig_du  = [];      % 4 x 1  acceptable correction per channel          -> R_du
    opts.sig_ddu = [];      % 4 x 1  acceptable rate of change of the correction -> S

    % IPOPT
    opts.max_iter    = 500;
    opts.tol         = 1e-8;
    opts.print_level = 5;   % 0 = silent, 5 = standard IPOPT log

end


function sg = deriveSigmas(mp, sig_p, tau_c, kappa_du)
% DERIVESIGMAS  Bryson tolerances from ONE position tolerance and ONE time
%               scale.
%
%   sg.x   : 12 x 1 state tolerances [p_I; v_I; angles; rates]
%   sg.u   :  4 x 1 input tolerances [T; L; M; N]
%   sg.du  :  4 x 1 tolerated input correction
%   sg.ddu :  4 x 1 tolerated rate of change of the correction

    s_v  = sig_p / tau_c;                        % velocity  [m/s]
    s_th = s_v / (mp.g * tau_c);                 % tilt      [rad]
    s_w  = s_th / tau_c;                         % body rate [rad/s]

    s_T = mp.m   * s_v / tau_c;                  % collective thrust -> velocity error  [N]
    s_L = mp.Jxx * s_w / tau_c;                  % roll  torque      -> rate error      [N*m]
    s_M = mp.Jyy * s_w / tau_c;                  % pitch torque      -> rate error      [N*m]
    s_N = mp.Jzz * s_w / tau_c;                  % yaw   torque      -> rate error      [N*m]

    sg.x   = [sig_p * ones(3,1); s_v * ones(3,1); s_th * ones(3,1); s_w * ones(3,1)];
    sg.u   = [s_T; s_L; s_M; s_N];
    sg.du  = kappa_du * sg.u;
    sg.ddu = sg.du / tau_c;

end


function Ainv = allocationInverse(mp)
% ALLOCATIONINVERSE  [T1;T2;T3;T4] = Ainv * [T;L;M;N]  (configuration C,
%                    same rows as C00_FlatnessMap).

    a = sqrt(2) / (4*mp.b);
    c = 1 / (4*mp.Beta);
    Ainv = [ 1/4,  a,  a, -c;
             1/4, -a,  a,  c;
             1/4, -a, -a, -c;
             1/4,  a, -a,  c ];

end


function opts = mergeOptions(opts, user)
% MERGEOPTIONS  Overwrites default fields with user-provided ones.
%               Unknown fields raise an error (catches typos).

    f = fieldnames(user);
    for i = 1:numel(f)
        if ~isfield(opts, f{i})
            error('A02_GeneratePWC_reference:UnknownOption', 'Unknown option "%s".', f{i});
        end
        opts.(f{i}) = user.(f{i});
    end

end


function F_int = buildIntervalFunction(model, nx, nu, M, h, q_diag)
% BUILDINTERVALFUNCTION  Symbolic one-interval map (CasADi SX):
%
%   [x_next, J_k] = F_int(x_k, u_k, XR_k)
%
%   x_k   : nx x 1  state at the start of the interval
%   u_k   : nu x 1  input, constant over the interval (ZOH)
%   XR_k  : nx x M  continuous reference at the M substep instants
%   x_next: nx x 1  state after Ts (M fixed RK4 steps of size h)
%   J_k   : 1  x 1  running tracking cost (left rectangle rule on substeps)

    import casadi.*

    xs  = SX.sym('x',  nx);
    us  = SX.sym('u',  nu);
    xrs = SX.sym('xr', nx, M);

    f = @(x, u) model.dynamics(0, x, u);

    xk = xs;
    Jk = 0;
    for j = 1:M
        e  = xk - xrs(:, j);
        Jk = Jk + h * sum1(q_diag .* e.^2);

        k1 = f(xk,              us);
        k2 = f(xk + (h/2) * k1, us);
        k3 = f(xk + (h/2) * k2, us);
        k4 = f(xk +  h    * k3, us);
        xk = xk + (h/6) * (k1 + 2*k2 + 2*k3 + k4);
    end

    F_int = Function('F_int', {xs, us, xrs}, {xk, Jk}, ...
                     {'x', 'u', 'xr'}, {'x_next', 'J'});

end


function u_bar = computeMeanInput_ode45(traj_obj, mp, t_now, Ts, quad_options)
% COMPUTEMEANINPUT_ODE45  Exact (to solver tolerance) time-average of
%                         u(t) = A02_FlatnessMap.map(...) over
%                         [t_now, t_now+Ts], via the auxiliary ODE
%                         dz/dt = u(t), z(t_now)=0.
%
%   u_bar = z(t_now+Ts) / Ts

    z0    = zeros(4, 1);
    tspan = [t_now, t_now + Ts];

    [~, z_sol] = ode45(@(t, z) evalU(t, traj_obj, mp), ...
                        tspan, z0, quad_options);

    u_bar = z_sol(end, :)' / Ts;

end


function u = evalU(t, traj_obj, mp)
% EVALU  Evaluates the flatness-consistent input u(t) = [T;L;M;N] at
%        a single time instant. Used as the (state-independent) vector
%        field of the auxiliary quadrature ODE dz/dt = u(t).

    [s, ds, dds, ddds, dddds] = traj_obj.get_flat_outputs(t);
    [~, u] = A02_FlatnessMap.map(mp, s, ds, dds, ddds, dddds);

end