% B01_QPBUILDER  Sparse QP matrix factory for the 16-state B01 model, OSQP format.
%
% Builds the standard OSQP form:
%
%   minimize   0.5 * z' * P * z + q' * z
%   subject to l <= A * z <= u
%
% Decision vector (SAME ordering convention as QPBuilder.m, nx=16, nu=4):
%   z = [x_0', u_0', x_1', u_1', ..., u_{N-1}', x_N']'
%
% Variables are ABSOLUTE (x_k, u_k), not deviations. 
% The linear cost term f = -2*Q*xref carries
% the reference-tracking behaviour, so (x-xref)'Q(x-xref) is recovered
% up to a constant (dropped, does not affect the argmin).
%
% Constraint stacking: equality rows (initial condition + LTV dynamics)
% on top, box rows (state/input bounds) below l=u on the equality block encodes 
% equality directly in OSQP's l<=Az<=u form.

classdef B01_QPBuilder
    methods (Static)

        function [P, q, A, l, u] = build(Ad_seq, Bd_seq, d_seq, ...
                                          x0_meas, xref_seq, uref_seq, ...
                                          Q, R, Qf, u_lims, x_lims)
            % Args: nx=16, nu=4 fixed for B01.
            % Returns: OSQP-format (P,q,A,l,u) instead of quadprog-format.

            N  = length(Ad_seq);
            nx = size(Ad_seq{1}, 1);   % 16
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

            % x_0 = x_meas
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

            %% BOX BLOCK (state/input bounds), OSQP style: l <= I*z <= u
            lb = -inf(nz, 1);
            ub =  inf(nz, 1);

            for k = 1:N
                idx_xk = (k-1)*(nx+nu) + 1;
                idx_uk = idx_xk + nx;

                if ~isempty(x_lims) && k > 1
                    lb(idx_xk : idx_xk+nx-1) = x_lims.min;
                    ub(idx_xk : idx_xk+nx-1) = x_lims.max;
                end

                lb(idx_uk : idx_uk+nu-1) = u_lims.min;
                ub(idx_uk : idx_uk+nu-1) = u_lims.max;
            end

            idx_xN = N*(nx+nu) + 1;
            if ~isempty(x_lims)
                lb(idx_xN : idx_xN+nx-1) = x_lims.min;
                ub(idx_xN : idx_xN+nx-1) = x_lims.max;
            end

            %% STACK INTO OSQP FORM
            % A = [Aeq; I],  l = [beq; lb],  u = [beq; ub]
            % Equality rows have l==u by construction (beq on both sides).
            A = [Aeq; speye(nz)];
            l = [beq; lb];
            u = [beq; ub];

        end
    end
end