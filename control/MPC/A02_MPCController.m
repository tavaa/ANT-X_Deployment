% A02_MPCCONTROLLER  LTV-MPC for the 12-state A02 model (deployment).
%
% Architecture:
%   1. Nominal input sequence u_hat: at the first call the reference input,
%      afterwards the previous optimal sequence shifted by
%      one step with the last element repeated "warm start".
%   2. Nominal state trajectory x_hat: nonlinear propagation (ode45) of
%      u_hat from the initial state of the prediction.
%   3. Linearization at every (x_hat_k, u_hat_k) and exact ZOH discretization.
%   4. Affine term d_k = x_hat_{k+1} - Ad_k*x_hat_k - Bd_k*u_hat_k:
%   5. QP (A02_QPBuilder, OSQP): tracking cost on (x - x_ref), (u - u_ref),
%      terminal weight Qf, state box, motor polytope on u = [T;L;M;N].
%   6. Only the first optimal input is returned (receding horizon).
%
% Weights: Bryson rule, same sigmas as the A02 NLP (A02_GeneratePWC_reference):
% or custom Q, Q_f, R -> from Control Parameters
% Constraints: state box = ControlParameters x_min_12 / x_max_12,
%              inputs     = physical rotor limits u_min_16 <= ainv*u <= u_max_16.

classdef A02_MPCController < handle

    properties
        model       % A02_SimplifiedModel_Deployment instance (12-state)
        cp          % ControlParameters

        Q           % [12x12] state error penalty
        R           % [4x4]   input error penalty
        Qf          % [12x12] terminal state penalty
        
        % optional
        sig_x       % [12x1]  state tolerances (Bryson)
        sig_u       % [4x1]   input tolerances (Bryson)

        Ainv        % [4x4]   [T;L;M;N] -> [T1..T4]

        nominal_x   % [12 x N+1] nominal state trajectory (linearization points)
        nominal_u   % [4  x N]   nominal input trajectory (warm start)

        ode_opts    % odeset for the nonlinear rollout / one-step prediction
        osqp_eps    % OSQP eps_abs = eps_rel

        osqp_prob         % persistent OSQP object (update-based warm start)
        osqp_initialized  % true after the first setup()
    end

    methods

        % CONSTRUCTOR
        %   model_obj : A02_SimplifiedModel_Deployment(mp)
        %   sig_p     : (optional) position tolerance [m], default 0.01 (same as NLP)
        %   tau_c     : (optional) time scale [s], default 2*Ts (same as NLP)
        function obj = A02_MPCController(model_obj, sig_p, tau_c)
            obj.model = model_obj;
            obj.cp    = ControlParameters;
            mp        = model_obj.mp;

            if nargin < 2 || isempty(sig_p), sig_p = 0.01;           end
            if nargin < 3 || isempty(tau_c), tau_c = 2*obj.cp.Ts;   end

            [obj.sig_x, obj.sig_u] = A02_MPCController.bryson_sigmas(mp, sig_p, tau_c);
            
            % A: use Bryson's sigma for Q, Q_f, R.
            %obj.Q  = diag(1 ./ obj.sig_x.^2);
            %obj.R  = diag(1 ./ obj.sig_u.^2);
            %obj.Qf = 10 * obj.Q;

            % B: custom weights from ControlParameters
            obj.Q  = diag(obj.cp.Q_diag_12(:));
            obj.R  = diag(obj.cp.R_diag_12(:));
            obj.Qf = 10 * obj.Q;
            
            % inverse: [T_1, ... , T_4]' = Ainv * [T, L, M, N]'
            a = sqrt(2) / (4*mp.b);
            c = 1 / (4*mp.Beta);
            obj.Ainv = [ 1/4,  a,  a, -c;
                         1/4, -a,  a,  c;
                         1/4, -a, -a, -c;
                         1/4,  a, -a,  c ];

            obj.ode_opts = odeset('RelTol', 1e-7, 'AbsTol', 1e-9);
            obj.osqp_eps = 1e-5;

            obj.osqp_initialized = false;
        end

        % INIT  Cold-start nominal trajectory from the (12-state) reference.
        %
        %   x0    - [12x1] initial state of the prediction
        %   x_ref - [12xT] reference states, starting at the instant of the
        %                  first prediction (k+1 in one-step-delay mode)
        %   u_ref - [4xT]  reference inputs, same alignment
        function init(obj, x0, x_ref, u_ref)
            N  = obj.cp.p;
            nx = 12;
            nu = 4;

            T_ref = size(x_ref, 2);

            obj.nominal_x = zeros(nx, N+1);
            obj.nominal_u = zeros(nu, N);
            
            % populate nominal trajectory with
            % reference trajectory (MPC init)
            for k = 1:N+1
                obj.nominal_x(:, k) = x_ref(:, min(k, T_ref));
            end
            for k = 1:N
                obj.nominal_u(:, k) = u_ref(:, min(k, size(u_ref, 2)));
            end

            obj.nominal_x(:, 1) = x0;
        end

        % SOLVE  Receding horizon optimization step.
        %
        %   x_meas      - [12x1]     current (measured) state
        %   xref_seq    - [12x(N+1)] state reference over the horizon, starting
        %                            at the instant of the prediction's first state
        %   uref_seq    - [4xN]      input reference over the horizon
        %   segments    - [5xm] (optional) [duration; u] of the commands that
        %                 act between the time stamp of x_meas and t_app, in
        %                 time order. The prediction starts at t_app from
        %                 x_hat = Phi_segments(x_meas) and the returned input is
        %                 meant for [t_app, t_app+Ts]. A [4x1] u is read as
        %                 [Ts; u]. If omitted/empty: prediction from x_meas.
        %
        % Returns:
        %   u_next     - [4x1] optimal input [T; L; M; N]
        %   debug_info - struct: x_opt, u_opt, x0_pred, status
        function [u_next, debug_info] = solve(obj, x_meas, xref_seq, uref_seq, segments)
            if nargin < 5, segments = []; end

            N  = obj.cp.p;
            Ts = obj.cp.Ts;
            nx = 12;
            nu = 4;

            %% Initial state of the prediction
            x0 = x_meas(:);
            %display(segments, 2) % debug
            if ~isempty(segments)
                if size(segments, 1) == 4            % old form: one input over Ts
                    segments = [Ts; segments(:, 1)];
                end
                for j = 1:size(segments, 2)
                    dur = segments(1, j);
                    if dur > 1e-9
                        u_j = segments(2:5, j);
                        [~, xs] = ode45(@(t, x) obj.model.dynamics(t, x, u_j), [0 dur], x0, obj.ode_opts);
                        x0 = xs(end, :)';
                    end
                end
            end

            %% Yaw branch of the reference aligned to the initial state
            xref_seq = A02_MPCController.align_yaw(xref_seq, x0(9));

            obj.nominal_x(:, 1) = x0;

            %% Sequential LTV linearization along the nominal trajectory
            Ad_seq = cell(1, N);
            Bd_seq = cell(1, N);
            d_seq  = cell(1, N);

            x_curr_sim = x0;

            for k = 1:N
                u_lin = obj.nominal_u(:, k);
                x_lin = x_curr_sim;

                [Ac, Bc] = obj.model.get_jacobians(x_lin, u_lin);
                [Ad, Bd] = Discretization.discretize(Ac, Bc, Ts, 'zoh');   % full-system exact ZOH

                % Nonlinear propagation for the affine correction term
                [~, x_sim_temp] = ode45(@(t, x) obj.model.dynamics(t, x, u_lin), [0 Ts], x_lin, obj.ode_opts);
                x_next_nl = x_sim_temp(end, :)';

                d_val = x_next_nl - (Ad * x_lin + Bd * u_lin);

                % Fixed, fully-dense sparsity pattern on Ad, Bd (OSQP 'Ax' update)
                EPS_PATTERN = 1e-12;
                Ad = Ad + EPS_PATTERN;
                Bd = Bd + EPS_PATTERN;

                Ad_seq{k} = Ad;
                Bd_seq{k} = Bd;
                d_seq{k}  = d_val;

                x_curr_sim = x_next_nl;
                obj.nominal_x(:, k+1) = x_next_nl;
            end

            %% Constraints
            u_lims.A   = obj.Ainv;              % rotor thrusts from [T;L;M;N]
            u_lims.min = obj.cp.u_min_16;       % physical rotor limits [N]
            u_lims.max = obj.cp.u_max_16;
            x_lims.min = obj.cp.x_min_12;
            x_lims.max = obj.cp.x_max_12;

            %% Build OSQP-format QP
            [P, q, A, l, u] = A02_QPBuilder.build( ...
                Ad_seq, Bd_seq, d_seq, ...
                x0, xref_seq, uref_seq, ...
                obj.Q, obj.R, obj.Qf, u_lims, x_lims);

            %% OSQP solve, setup once, then update
            if ~obj.osqp_initialized
                obj.osqp_prob = osqp;
                obj.osqp_prob.setup(P, q, A, l, u, ...
                    'warm_start', true, 'verbose', false, 'polish', false, ...
                    'eps_abs', obj.osqp_eps, 'eps_rel', obj.osqp_eps);
                obj.osqp_initialized = true;
            else
                Ax_new = nonzeros(A);
                obj.osqp_prob.update('q', q, 'l', l, 'u', u, 'Ax', Ax_new);
            end

            res = obj.osqp_prob.solve();
            debug_info.x0_pred = x0;
            debug_info.status  = res.info.status;

            if res.info.status_val ~= 1
                warning('MPC:QP_Fail', 'OSQP solver failure (status: %s). Reverting to the nominal (warm-start) input.', ...
                    res.info.status);
                u_next = obj.nominal_u(:, 1);

                obj.nominal_u = [obj.nominal_u(:, 2:end), obj.nominal_u(:, end)];

                debug_info.x_opt = obj.nominal_x;
                debug_info.u_opt = obj.nominal_u;
                return;
            end

            z_opt = res.x;

            %% De-stacking
            x_opt_seq = zeros(nx, N+1);
            u_opt_seq = zeros(nu, N);

            for k = 0:N-1
                idx_x = k*(nx+nu) + 1;
                idx_u = idx_x + nx;

                x_opt_seq(:, k+1) = z_opt(idx_x : idx_x+nx-1);
                u_opt_seq(:, k+1) = z_opt(idx_u : idx_u+nu-1);
            end

            idx_xN = N*(nx+nu) + 1;
            x_opt_seq(:, N+1) = z_opt(idx_xN : idx_xN+nx-1);

            u_next = u_opt_seq(:, 1);

            %% Warm-start shift (last element repeated)
            obj.nominal_u = [u_opt_seq(:, 2:end), u_opt_seq(:, end)];
            obj.nominal_x = [x_opt_seq(:, 2:end), x_opt_seq(:, end)];

            debug_info.x_opt = x_opt_seq;
            debug_info.u_opt = u_opt_seq;
        end

        % NEXT_PLANNED_INPUT  Input the current plan intends for the interval that
        %                      starts at the next t_app (warm-start element 1).
        %                      Useful as a fallback BEFORE calling solve().
        function u = next_planned_input(obj)
            u = obj.nominal_u(:, 1);
        end

    end

    methods (Static)

        function [sig_x, sig_u] = bryson_sigmas(mp, sig_p, tau_c)
            % BRYSON_SIGMAS  Same rule as A02_GeneratePWC_reference (deriveSigmas).
            s_v  = sig_p / tau_c;                % velocity  [m/s]
            s_th = s_v / (mp.g * tau_c);         % tilt      [rad]
            s_w  = s_th / tau_c;                 % body rate [rad/s]

            sig_x = [sig_p*ones(3,1); s_v*ones(3,1); s_th*ones(3,1); s_w*ones(3,1)];
            sig_u = [mp.m*s_v; mp.Jxx*s_w; mp.Jyy*s_w; mp.Jzz*s_w] / tau_c;
        end

    end

    methods (Static, Access = private)

        function xr = align_yaw(xr, psi0)
            % ALIGN_YAW  Shift the reference yaw row by 2*pi*n so that its
            %            first sample is within +-pi of psi0 (continuity of
            %            the window is preserved).
            n = round((psi0 - xr(9, 1)) / (2*pi));
            xr(9, :) = xr(9, :) + 2*pi*n;
        end

    end
end
