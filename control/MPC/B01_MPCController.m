% B01_MPCCONTROLLER  LTV-MPC for the 16-state B01 model.
%
% Implementation notes:
%
%   Reference padding: the reference generator produces a 12-state
%   (no-delay) trajectory. This controller pads it to 16 states
%   internally using u_ref for the T1..T4 slots (T_i,ref := T_i,ref^d),
%   since the reference has no motor-delay states to draw from
%   directly.
%
%   Discretization: full exact ZOH (matrix exponential, via
%   Discretization.discretize(...,'zoh')) on the entire 16-state/
%   4-input system at every linearization step.
%
%   OSQP as the solver for the quadratic optimization problem

classdef B01_MPCController < handle

    properties
        model       % B01_ANTX_quadcopter instance, delay_motors=true (16-state)
        cp          % ControlParameters

        Q           % [16x16] State error penalty matrix
        R           % [4x4]   Input penalty matrix
        Qf          % [16x16] Terminal state penalty matrix

        nominal_x   % [16 x N+1] Nominal state trajectory for linearization
        nominal_u   % [4  x N]   Nominal input trajectory for linearization

        osqp_prob         % persistent OSQP problem object (update-based warm start)
        osqp_initialized  % true after the first prob.setup() call
    end

    methods

        % CONSTRUCTOR
        function obj = B01_MPCController(model_obj)
            % model_obj : B01_ANTX_quadcopter(mp, true) -- 16-state, delayed
            obj.model = model_obj;
            obj.cp    = ControlParameters;

            obj.Q  = diag(obj.cp.Q_diag_16);
            obj.Qf = diag(obj.cp.Qf_diag_16);
            obj.R  = diag(obj.cp.R_diag_16);

            obj.osqp_initialized = false;
        end

        % INIT  Cold-start nominal trajectory from a 12-state (no-delay)
        % reference, padded to 16 states.
        %
        %   x0    - [16x1] Initial measured state
        %   x_ref - [12xT] Reference state trajectory
        %   u_ref - [4xT]  Reference input trajectory
        function init(obj, x0, x_ref, u_ref)
            N  = obj.cp.p;
            nx = 16;
            nu = 4;

            x_ref_16 = B01_MPCController.pad_reference(x_ref, u_ref);
            T_ref = size(x_ref_16, 2);

            obj.nominal_x = zeros(nx, N+1);
            obj.nominal_u = zeros(nu, N);

            for k = 1:N+1
                idx = min(k, T_ref);
                obj.nominal_x(:, k) = x_ref_16(:, idx);
            end

            for k = 1:N
                idx = min(k, size(u_ref, 2));
                obj.nominal_u(:, k) = u_ref(:, idx);
            end

            obj.nominal_x(:, 1) = x0;
        end

        % SOLVE  Receding horizon optimization step.
        %
        %   x_meas   - [16x1]     Current state feedback measurement
        %   xref_seq - [12x(N+1)] Reference state trajectory over horizon (no-delay, 12-state)
        %   uref_seq - [4xN]      Reference input trajectory over horizon
        %
        % Returns:
        %   u_next     - [4x1]  Optimal desired motor thrusts [T1_d;...;T4_d]
        %   debug_info
        function [u_next, debug_info] = solve(obj, x_meas, xref_seq, uref_seq)
            N  = obj.cp.p;
            Ts = obj.cp.Ts;
            nx = 16;
            nu = 4;

            xref_seq_16 = B01_MPCController.pad_reference(xref_seq, uref_seq);

            obj.nominal_x(:, 1) = x_meas;

            % Sequential LTV linearization loop
            Ad_seq = cell(1, N);
            Bd_seq = cell(1, N);
            d_seq  = cell(1, N);

            x_curr_sim = x_meas;

            for k = 1:N
                u_lin = obj.nominal_u(:, k);
                x_lin = x_curr_sim;

                [Ac, Bc] = obj.model.get_jacobians(x_lin, u_lin);
                [Ad, Bd] = Discretization.discretize(Ac, Bc, Ts, 'zoh'); % Full-system exact ZOH.

                % Nonlinear propagation for the affine correction term
                [~, x_sim_temp] = ode45(@(t,x) obj.model.dynamics(t, x, u_lin), [0 Ts], x_lin);
                x_next_nl = x_sim_temp(end, :)';

                d_val = x_next_nl - (Ad * x_lin + Bd * u_lin);

                % Force a FIXED, fully-dense sparsity pattern on Ad, Bd
                % before they reach B01_QPBuilder
                EPS_PATTERN = 1e-12;
                Ad = Ad + EPS_PATTERN;
                Bd = Bd + EPS_PATTERN;

                Ad_seq{k} = Ad;
                Bd_seq{k} = Bd;
                d_seq{k}  = d_val;

                x_curr_sim = x_next_nl;
                obj.nominal_x(:, k+1) = x_next_nl;
            end

            % Constraints
            u_lims.min = obj.cp.u_min_16;
            u_lims.max = obj.cp.u_max_16;
            x_lims.min = obj.cp.x_min_16;
            x_lims.max = obj.cp.x_max_16;

            % Build OSQP-format QP
            [P, q, A, l, u] = B01_QPBuilder.build(...
                Ad_seq, Bd_seq, d_seq, ...
                x_meas, xref_seq_16, uref_seq, ...
                obj.Q, obj.R, obj.Qf, u_lims, x_lims);

            % OSQP solve -- setup once.
            if ~obj.osqp_initialized
                obj.osqp_prob = osqp;
                obj.osqp_prob.setup(P, q, A, l, u, ...
                    'warm_start', true, 'verbose', false, 'polish', false);
                obj.osqp_initialized = true;
            else
                Ax_new = nonzeros(A);
                obj.osqp_prob.update('q', q, 'l', l, 'u', u, 'Ax', Ax_new);
            end

            res = obj.osqp_prob.solve();

            if res.info.status_val ~= 1
                warning('MPC:QP_Fail', 'OSQP solver failure (status: %s). Reverting to feedforward input.', ...
                    res.info.status);
                u_next = obj.nominal_u(:, 1);

                obj.nominal_u = [obj.nominal_u(:, 2:end), obj.nominal_u(:, end)];

                debug_info.x_opt = obj.nominal_x;
                debug_info.u_opt = obj.nominal_u;
                return;
            end

            z_opt = res.x;

            % De-stacking 
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

            % Warm-start shift
            obj.nominal_u = [u_opt_seq(:, 2:end), u_opt_seq(:, end)];
            obj.nominal_x = [x_opt_seq(:, 2:end), x_opt_seq(:, end)];

            debug_info.x_opt = x_opt_seq;
            debug_info.u_opt = u_opt_seq;
        end

    end

    methods (Static, Access = private)

        function x16 = pad_reference(x12, uref)
            % Pads a 12-state (no-delay) reference to 16 states by
            % appending T1..T4 = corresponding uref column (or the last
            % available uref column, if x12 has one more column than
            % uref).
            M  = size(x12, 2);
            Mu = size(uref, 2);

            x16 = zeros(16, M);
            x16(1:12, :) = x12;

            for k = 1:M
                idx = min(k, Mu);
                x16(13:16, k) = uref(:, idx);
            end
        end

    end
end