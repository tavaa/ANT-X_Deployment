% A02_QPBUILDER  Sparse QP matrix factory for the 12-state A02 model, OSQP format.
%
% Builds the standard OSQP form:
%
%   minimize   0.5 * z' * P * z + q' * z
%   subject to l <= A * z <= u
%
% Decision vector (nx=12, nu=4):
%   z = [x_0', u_0', x_1', u_1', ..., u_{N-1}', x_N']'
%   u_k = [T; L; M; N]
%
% Variables are ABSOLUTE (x_k, u_k).
% The linear cost term f = -2*Q*xref carries the reference-tracking
% behaviour, so (x-xref)'Q(x-xref) is recovered up to a constant.
%
% Constraint stacking: equality rows (initial condition + LTV dynamics)
% on top, box rows (state bounds) below, motor-polytope rows at the bottom:
%   T_min <= Ainv * u_k <= T_max      (single-rotor limits, linear in u_k)
% l = u on the equality block encodes equality directly in OSQP's l <= Az <= u form.

classdef A02_QPBuilder
    methods (Static)

        function [P, q, A, l, u] = build(Ad_seq, Bd_seq, d_seq, ...
                                          x0_meas, xref_seq, uref_seq, ...
                                          Q, R, Qf, u_lims, x_lims)
            % Args (nx = 12, nu = 4 for A02):
            %   Ad_seq, Bd_seq, d_seq : 1xN cells, LTV model x+ = Ad x + Bd u + d
            %   x0_meas               : nx x 1   initial state of the prediction
            %   xref_seq              : nx x (N+1) state reference over the horizon
            %   uref_seq              : nu x N     input reference over the horizon
            %   Q, R, Qf              : weights
            %   u_lims                : struct, fields A (4x4 = Ainv), min, max (4x1 rotor limits [N])
            %   x_lims                : struct with fields min/max
            % Returns: OSQP-format (P, q, A, l, u).

            N  = length(Ad_seq);
            nx = size(Ad_seq{1}, 1);   % 12
            nu = size(Bd_seq{1}, 2);   % 4
            nz = (N+1)*nx + N*nu;

            %% COST (P, q)
            H_blocks = cell(1, 2*N + 1);
            f_cell   = cell(1, 2*N + 1);

            Q_sp  = sparse(Q);
            R_sp  = sparse(R);
            Qf_sp = sparse(Qf);

            for k = 1:N
                H_blocks{2*k - 1} = 2 * Q_sp;
                H_blocks{2*k}     = 2 * R_sp;

                f_cell{2*k - 1} = -2 * Q * xref_seq(:, k);
                f_cell{2*k}     = -2 * R * uref_seq(:, k);
            end

            H_blocks{2*N + 1} = 2 * Qf_sp;
            f_cell{2*N + 1}   = -2 * Qf * xref_seq(:, N+1);

            P = blkdiag(H_blocks{:});
            q = vertcat(f_cell{:});

            %% EQUALITY BLOCK (initial condition + LTV dynamics)
            num_eq = (N+1) * nx;
            nz_est = num_eq * (nx + nu + 1);
            Aeq = spalloc(num_eq, nz, nz_est);
            beq = zeros(num_eq, 1);

            % x_0 = x_meas (or the predicted state at the application time)
            Aeq(1:nx, 1:nx) = speye(nx);
            beq(1:nx)       = x0_meas;

            row_idx = nx + 1;
            for k = 1:N
                idx_xk   = (k-1)*(nx+nu) + 1;
                idx_uk   = idx_xk + nx;
                idx_xkp1 = idx_uk + nu;

                rows = row_idx : row_idx+nx-1;

                Aeq(rows, idx_xk   : idx_xk+nx-1)   = sparse(-Ad_seq{k});
                Aeq(rows, idx_uk   : idx_uk+nu-1)   = sparse(-Bd_seq{k});
                Aeq(rows, idx_xkp1 : idx_xkp1+nx-1) = speye(nx);

                beq(rows) = d_seq{k};

                row_idx = row_idx + nx;
            end

            %% BOX BLOCK (state bounds), OSQP style: l <= I*z <= u
            % Inputs have no box: their physical limit is the polytope below.
            lb = -inf(nz, 1);
            ub =  inf(nz, 1);

            for k = 1:N
                idx_xk = (k-1)*(nx+nu) + 1;

                if ~isempty(x_lims) && k > 1          % x_0 is fixed, not bounded
                    lb(idx_xk : idx_xk+nx-1) = x_lims.min;
                    ub(idx_xk : idx_xk+nx-1) = x_lims.max;
                end
            end

            idx_xN = N*(nx+nu) + 1;
            if ~isempty(x_lims)
                lb(idx_xN : idx_xN+nx-1) = x_lims.min;
                ub(idx_xN : idx_xN+nx-1) = x_lims.max;
            end

            %% MOTOR POLYTOPE BLOCK: T_min <= Ainv * u_k <= T_max
            nc   = size(u_lims.A, 1);                 % 4 rotors
            Apol = spalloc(N*nc, nz, N*nc*nu);
            lpol = zeros(N*nc, 1);
            upol = zeros(N*nc, 1);

            for k = 1:N
                idx_uk = (k-1)*(nx+nu) + nx + 1;
                rows   = (k-1)*nc + (1:nc);

                Apol(rows, idx_uk : idx_uk+nu-1) = sparse(u_lims.A);
                lpol(rows) = u_lims.min;
                upol(rows) = u_lims.max;
            end

            %% STACK INTO OSQP FORM
            % A = [Aeq; I; Apol],  l = [beq; lb; lpol],  u = [beq; ub; upol]
            % Equality rows have l==u by construction (beq on both sides).
            A = [Aeq; speye(nz); Apol];
            l = [beq; lb; lpol];
            u = [beq; ub; upol];

        end
    end
end
