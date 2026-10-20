classdef B01_FlatnessMap
    % B01_FlatnessMap  Differential flatness inverse map, 16-state model.
    %
    % Flat outputs: sigma = [x_I, y_I, z_I, psi] in R^4. Reconstructs the
    % reference state x in R^16 and input u in R^4 from the flat outputs
    % and their time derivatives, without integration.
    %
    % NOTE: this map does NOT include motor delay. The 4 motor-thrust
    % states are padded algebraically: T_i,ref := T_i^d,ref, i.e. the
    % reference assumes instantaneous actuation.
    %
    % Flat output smoothness requirements
    %   Position sigma(1:3) : C^4 (snap required for torque reconstruction)
    %   Yaw      sigma(4)   : C^2 (yaw acceleration required)
    %
    % Reference frames
    %   Inertial : NED  (North-East-Down)   e_z = [0;0;1] points Down
    %   Body     : FRD  (Forward-Right-Down)
    %   Thrust F_thrust^I = R * [0;0;-T]  (thrust along -z_B, positive upward)
    %
    % State reconstruction x = h(sigma, dsigma, ddsigma, dddsigma)
    %   x(1:3)   = p_I         Position directly from flat outputs.
    %   x(4:6)   = v_I         Inertial (NED) velocity, directly from dsigma.
    %   x(7:8)   = phi, theta  Derived from thrust vector (z_B) and desired heading.
    %   x(9)     = psi         Euler yaw extracted from R (psi_check).
    %   x(10:12) = p, q, r     Body rates derived from position jerk and yaw rate.
    %   x(13:16) = T1..T4      Individual motor thrusts (no-delay reference).
    %
    % Input reconstruction u = i(sigma, ..., ddddsigma)
    %   u(1:4) = T1..T4   desired motor thrusts (equal to x(13:16) here)
    %
    % Singularities
    %   T = 0 (free fall): thrust direction undefined.
    %   z_B parallel to x_C: gimbal-like singularity at 90 deg pitch aligned with yaw.

    methods (Static)

        function [xref, uref] = map(mp, sigma, dsigma, ddsigma, dddsigma, ddddsigma)
            % MAP  Algebraic reconstruction of reference state and input.
            %
            %   [xref, uref] = B01_FlatnessMap.map(mp, sigma, dsigma, ddsigma, dddsigma, ddddsigma)
            %
            %   Inputs
            %     mp         : ModelParameters (fields m,g,Jxx,Jyy,Jzz,b,Beta,Lp,Mq,Nr)
            %     sigma      : [4x1] flat outputs     [x_I; y_I; z_I; psi]
            %     dsigma     : [4x1] 1st derivatives  [vx_I; vy_I; vz_I; psi_dot]
            %     ddsigma    : [4x1] 2nd derivatives  [ax_I; ay_I; az_I; psi_ddot]
            %     dddsigma   : [4x1] 3rd derivatives  [jx_I; jy_I; jz_I; psi_dddot] (jerk)
            %     ddddsigma  : [4x1] 4th derivatives  [sx_I; sy_I; sz_I; psi_ddddot] (snap)
            %
            %   Outputs
            %     xref : [16x1] reference state  [p_I; v_I; phi; theta; psi; p; q; r; T1; T2; T3; T4]
            %     uref : [4x1]  reference input   [T1; T2; T3; T4]

            %% Parameters
            g  = mp.g;
            m  = mp.m;

            Jxx = mp.Jxx;
            Jyy = mp.Jyy;
            Jzz = mp.Jzz;
            J = diag([Jxx Jyy Jzz]);

            b    = mp.b;
            beta = mp.Beta;
            Lp   = mp.Lp;
            Mq   = mp.Mq;
            Nr   = mp.Nr;


            %% Flat Outputs
            psi     = sigma(4);
            psi_dot = dsigma(4);
            psi_dd  = ddsigma(4);

            p_I = sigma(1:3);
            v_I = dsigma(1:3);
            a_I = ddsigma(1:3);
            j_I = dddsigma(1:3);
            s_I = ddddsigma(1:3);


            %% Net force vector in NED frame
            %
            %   m*a_I = m*g*e_z + R*[0;0;-T]
            %   => F_net := m*(a_I - [0;0;g]) = -T*z_B^I
            %
            % g appears with + sign because NED convention: e_z=[0;0;1]
            % points Down, and gravity IS m*g*e_z.

            F_net = m * (a_I - [0; 0; g]);
            T_ref = norm(F_net);

            assert(T_ref > 1e-4, ...
                'FlatnessMap:ZeroThrust', ...
                'T = %.4f N: trajectory passes through free fall (T=0). Increase thrust margin.', T_ref);

            z_B = -F_net / T_ref;


            %% Euler Angles [phi, theta, psi]
            %
            % psi is the 4th flat output (known exactly). phi, theta are
            % derived from z_B^I and the desired yaw:
            %   1. x_C = [cos(psi), sin(psi), 0]  (horizontal heading vector)
            %   2. y_B^I = (z_B^I x x_C) / ||(z_B^I x x_C)||
            %   3. x_B^I = y_B^I x z_B^I
            %   4. R = [x_B^I | y_B^I | z_B^I]
            %   5. Extract phi, theta, psi_check from R (Z-Y-X decomposition)

            x_C = [cos(psi); sin(psi); 0];
            y_B_cross = cross(z_B, x_C);
            y_B_norm  = norm(y_B_cross);

            assert(y_B_norm > 1e-6, ...
                'FlatnessMap:AttitudeSingularity', ...
                'z_B parallel to x_C: pitch singularity near 90 deg aligned with yaw.');

            y_B = y_B_cross / y_B_norm;
            x_B = cross(y_B, z_B);

            R = [x_B, y_B, z_B];

            % Extract Euler angles (Z-Y-X convention)
            %   R(3,1) = -sin(theta)
            %   R(3,2) =  sin(phi)*cos(theta)
            %   R(3,3) =  cos(phi)*cos(theta)
            psi_check = atan2(R(2,1), R(1,1));
            theta = asin(-R(3,1));
            phi   = atan2(R(3,2), R(3,3));

            yaw_err = atan2( ...
                sin(psi - psi_check), ...
                cos(psi - psi_check));
            %#ok<NASGU> kept for future diagnostic use, not asserted

            cthe = cos(theta);

            assert(abs(cthe) > 1e-6, ...
                'FlatnessMap:GimbalLock', ...
                'cos(theta) ~ 0.');


            %% Inertial velocity -- direct, no rotation needed
            %
            % Velocity is stored natively in NED, so v_I from dsigma IS
            % the state -- no rotation required.
            v_ref = v_I;


            %% Angular velocity and acceleration [p,q,r], [p_dot,q_dot,r_dot]
            %
            %   omega_hat     = R^T * R_dot
            %   omega_hat_dot = R^T * R_ddot - omega_hat^2

            F_dot = m * j_I;
            T_dot = -(z_B' * F_dot);

            z_B_dot = -(F_dot + T_dot*z_B) / T_ref;
            x_C_dot = psi_dot * [-sin(psi); cos(psi); 0];

            c_vec = cross(z_B, x_C);
            c_dot = cross(z_B_dot, x_C) + cross(z_B, x_C_dot);

            y_B_dot = (eye(3) - y_B*y_B') * c_dot / norm(c_vec);
            x_B_dot = cross(y_B_dot, z_B) + cross(y_B, z_B_dot);

            R_dot = [x_B_dot y_B_dot z_B_dot];
            omega_hat = R' * R_dot;

            p_ref = omega_hat(3,2);
            q_ref = omega_hat(1,3);
            r_ref = omega_hat(2,1);

            omega = [p_ref; q_ref; r_ref];

            F_ddot = m * s_I;
            T_ddot = -(z_B_dot' * F_dot + z_B' * F_ddot);

            z_B_ddot = -(F_ddot + T_ddot*z_B + 2*T_dot*z_B_dot) / T_ref;
            x_C_ddot = psi_dd * [-sin(psi); cos(psi); 0] + psi_dot^2 * [-cos(psi); -sin(psi); 0];
            c_ddot = cross(z_B_ddot, x_C) + 2*cross(z_B_dot, x_C_dot) + cross(z_B, x_C_ddot);

            P_y = eye(3) - y_B*y_B';

            y_B_ddot = P_y * c_ddot / norm(c_vec) - ( ...
                2 * y_B_dot * (y_B' * c_dot) ...
                + y_B * (y_B_dot' * c_dot) ...
                ) / norm(c_vec);

            x_B_ddot = ...
                cross(y_B_ddot, z_B) ...
                + 2*cross(y_B_dot, z_B_dot) ...
                + cross(y_B, z_B_ddot);

            R_ddot = [x_B_ddot y_B_ddot z_B_ddot];

            omega_hat_dot = R' * R_ddot - omega_hat * omega_hat;

            p_dot_ref = omega_hat_dot(3,2);
            q_dot_ref = omega_hat_dot(1,3);
            r_dot_ref = omega_hat_dot(2,1);

            omega_dot = [p_dot_ref; q_dot_ref; r_dot_ref];


            %% Torques [L, M, N] -- rigid body part (no damping yet)
            %
            %   tau = J*omega_dot + omega x (J*omega)
            %
            % tau is the moment WITHOUT the damping terms.

            tau = J*omega_dot + cross(omega, J*omega);


            %% Damping correction
            %
            % Dynamics: Jxx*p_dot = (Jyy-Jzz)*q*r + L_virt + Lp*p, so the
            % VIRTUAL moment (motor contribution only) is:
            %   L_virt = tau(1) - Lp*p
            %   M_virt = tau(2) - Mq*q
            %   N_virt = tau(3) - Nr*r
            % p,q,r already available (omega above) -- purely algebraic.

            L_virt = tau(1) - Lp*p_ref;
            M_virt = tau(2) - Mq*q_ref;
            N_virt = tau(3) - Nr*r_ref;


            %% Invert allocation -> individual motor thrusts
            %
            % Forward map:
            %   T = T1+T2+T3+T4
            %   L = (b/sqrt2)*( T1+T2-T3-T4)
            %   M = (b/sqrt2)*( T1-T2-T3+T4)
            %   N = beta     *( T1-T2+T3-T4)
            %
            % Rows are mutually orthogonal -> closed-form inverse, no
            % numerical matrix inversion needed.

            T1 = T_ref/4 + (sqrt(2)/(4*b))*(L_virt+M_virt) + N_virt/(4*beta);
            T2 = T_ref/4 + (sqrt(2)/(4*b))*(L_virt-M_virt) - N_virt/(4*beta);
            T3 = T_ref/4 - (sqrt(2)/(4*b))*(L_virt+M_virt) + N_virt/(4*beta);
            T4 = T_ref/4 - (sqrt(2)/(4*b))*(L_virt-M_virt) - N_virt/(4*beta);


            %% ASSEMBLE STATE AND INPUT REFERENCE VECTORS

            xref = [p_I;                    % 1-3   position NED
                    v_ref;                  % 4-6   inertial velocity
                    phi; theta; psi_check;  % 7-9   Euler angles
                    p_ref; q_ref; r_ref;    % 10-12 body rates
                    T1; T2; T3; T4];        % 13-16 motor thrusts (no-delay reference)

            % No motor delay in this map -> T_i,ref = T_i^d,ref.
            uref = [T1; T2; T3; T4];

        end
    end
end