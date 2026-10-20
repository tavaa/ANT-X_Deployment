classdef A01_FlatnessMap
% FLATNESSMAP  Differential flatness inverse map for the 13-state quaternion model.
% Flatness map for A01_SimplifiedModel_quat.m.
%
% Euler angle extraction (phi, theta, psi) is replaced
% by a quaternion extraction from R using Shepperd's method.
%
% Flat outputs
%   sigma = [x_I; y_I; z_I; psi]  in R^4
%
% Flat output smoothness requirements
%   Position sigma(1:3) : C^4 (snap required for torque reconstruction)
%   Yaw      sigma(4)   : C^2 (yaw acceleration required)
%
% Reference frames
%   Inertial : NED  (North-East-Down)
%   Body     : FRD  (Forward-Right-Down)
%   Thrust   : F_thrust^I = R * [0;0;-T]  (along -z_B, positive upward)
%
% State reconstruction
%   x(1:3)   = p_I         position from flat outputs
%   x(4:6)   = v_b         R^T * v_I
%   x(7:10)  = [q0;q1;q2;q3]  quaternion from R via Shepperd's method
%   x(11:13) = [p; q_r; r]    body rates from position jerk and yaw rate
%
% Input reconstruction
%   u(1) = T   from ||F_net||
%   u(2) = L   from snap + Euler equation
%   u(3) = M   from snap + Euler equation
%   u(4) = N   from snap + Euler equation
%
% Singularities
%   T = 0 (free fall): thrust direction undefined.
%   z_B parallel to x_C: same geometric singularity as A00.
%   No gimbal lock (quaternion representation is singularity-free).
 
    methods (Static)
 
        function [xref, uref] = map(mp, sigma, dsigma, ddsigma, dddsigma, ddddsigma)
        % MAP  Reconstruct reference state (13x1) and input (4x1).
        %
        %   Inputs and all internal computations are identical to A00_FlatnessMap.
        %   Only the attitude output is different: [phi;theta;psi] → [q0;q1;q2;q3].
        %
        %   Inputs
        %     mp        : ModelParameters
        %     sigma     : [4x1]  [x_I; y_I; z_I; psi]
        %     dsigma    : [4x1]  1st derivatives
        %     ddsigma   : [4x1]  2nd derivatives (acceleration)
        %     dddsigma  : [4x1]  3rd derivatives (jerk)
        %     ddddsigma : [4x1]  4th derivatives (snap)
        %
        %   Outputs
        %     xref : [13x1]  [p_I; v_b; q0; q1; q2; q3; p; q_r; r]
        %     uref : [4x1]   [T; L; M; N]
 
            %% Parameters
            g   = mp.g;    m   = mp.m;
            Jxx = mp.Jxx;  Jyy = mp.Jyy;  Jzz = mp.Jzz;
            J   = diag([Jxx Jyy Jzz]);
 

            %% Flat outputs
            psi     = sigma(4);
            psi_dot = dsigma(4);
            psi_dd  = ddsigma(4);
 
            p_I = sigma(1:3);
            v_I = dsigma(1:3);
            a_I = ddsigma(1:3);
            j_I = dddsigma(1:3);
            s_I = ddddsigma(1:3);
 

            %% Net force in NED frame
            %
            %   m*a_I = m*g*e3 + R*[0;0;-T]
            %   => F_net = m*(a_I - [0;0;g]) = -T*z_B^I
 
            F_net = m * (a_I - [0; 0; g]);
            T_ref = norm(F_net);
 
            assert(T_ref > 1e-4, ...
                'A01_FlatnessMap:ZeroThrust', ...
                'T = %.4f N: trajectory passes through free fall.', T_ref);
 
            z_B = -F_net / T_ref;
 

            %% Attitude geometric construction
            %
            %   x_C = [cos(psi); sin(psi); 0]
            %   y_B = (z_B x x_C) / ||z_B x x_C||
            %   x_B = y_B x z_B
            %   R   = [x_B | y_B | z_B]
 
            x_C       = [cos(psi); sin(psi); 0];
            y_B_cross = cross(z_B, x_C);
            y_B_norm  = norm(y_B_cross);
 
            assert(y_B_norm > 1e-6, ...
                'A01_FlatnessMap:AttitudeSingularity', ...
                'z_B parallel to x_C: pitch singularity near 90 deg aligned with yaw.');
 
            y_B = y_B_cross / y_B_norm;
            x_B = cross(y_B, z_B);
            R   = [x_B, y_B, z_B];
 
            %% Quaternion extraction from R (Shepperd's method)
            %
            % Numerically robust: chooses the formula with the largest
            % denominator to minimize round-off.
            % Convention: q0 >= 0.
 
            qq = RotationMatrix_to_quaternion(R);
            q0 = qq(1);  q1 = qq(2);  q2 = qq(3);  q3 = qq(4);
 
            %% Body frame velocities  v_b = R^T * v_I
            v_b = R' * v_I;
 
            %% Body angular rates from R_dot
            %
            %   omega_hat = R^T * R_dot
            %   p = omega_hat(3,2),  q_r = omega_hat(1,3),  r = omega_hat(2,1)
            %
            % R_dot is derived by differentiating z_B, y_B, x_B:
 
            F_dot   = m * j_I;
            T_dot   = -(z_B' * F_dot);
            z_B_dot = -(F_dot + T_dot*z_B) / T_ref;
            x_C_dot =  psi_dot * [-sin(psi); cos(psi); 0];
 
            c_vec   = cross(z_B, x_C);
            c_dot   = cross(z_B_dot, x_C) + cross(z_B, x_C_dot);
            y_B_dot = (eye(3) - y_B*y_B') * c_dot / norm(c_vec);
            x_B_dot = cross(y_B_dot, z_B) + cross(y_B, z_B_dot);
 
            R_dot    = [x_B_dot, y_B_dot, z_B_dot];
            omega_hat = R' * R_dot;
 
            p_ref  = omega_hat(3,2);
            q_r_ref= omega_hat(1,3);
            r_ref  = omega_hat(2,1);
            omega  = [p_ref; q_r_ref; r_ref];
 
            %% Angular accelerations from R_ddot
            %
            %   omega_hat_dot = R^T * R_ddot - omega_hat^2
            %   p_dot  = omega_hat_dot(3,2)
            %   q_r_dot= omega_hat_dot(1,3)
            %   r_dot  = omega_hat_dot(2,1)
 
            F_ddot   = m * s_I;
            T_ddot   = -(z_B_dot'*F_dot + z_B'*F_ddot);
            z_B_ddot = -(F_ddot + T_ddot*z_B + 2*T_dot*z_B_dot) / T_ref;
            x_C_ddot = psi_dd * [-sin(psi); cos(psi); 0] + psi_dot^2 * [-cos(psi); -sin(psi); 0];
            c_ddot   = cross(z_B_ddot,x_C) + 2*cross(z_B_dot,x_C_dot) + cross(z_B,x_C_ddot);
 
            P_y = eye(3) - y_B*y_B';
            y_B_ddot = P_y * c_ddot / norm(c_vec) ...
                     - (2*y_B_dot*(y_B'*c_dot) + y_B*(y_B_dot'*c_dot)) / norm(c_vec);
            x_B_ddot = cross(y_B_ddot,z_B) + 2*cross(y_B_dot,z_B_dot) + cross(y_B,z_B_ddot);
 
            R_ddot        = [x_B_ddot, y_B_ddot, z_B_ddot];
            omega_hat_dot = R' * R_ddot - omega_hat * omega_hat;
 
            omega_dot = [omega_hat_dot(3,2); omega_hat_dot(1,3); omega_hat_dot(2,1)];
 
            %% Torques  tau = J*omega_dot + omega x (J*omega)
            tau   = J * omega_dot + cross(omega, J*omega);
            L_ref = tau(1);
            M_ref = tau(2);
            N_ref = tau(3);
 
            %% Assemble state and input
            %
            %   State: [p_I; v_b; q0; q1; q2; q3; p; q_r; r]
 
            xref = [p_I;
                    v_b;
                    q0; q1; q2; q3;
                    p_ref; q_r_ref; r_ref];
 
            uref = [T_ref; L_ref; M_ref; N_ref];
        end
 
    end
end