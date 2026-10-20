classdef A00_FlatnessMap  
% FLATNESSMAP  Differential flatness inverse map for the 12-state body model.
% Flatness Map specifically built for models/A00_SimplifiedModel.m
% The quadrotor model is differentially
% flat with flat outputs:
%
% sigma = [x_I, y_I, z_I, psi]  in R^4
% This class reconstructs the full state x in R^12 and the input u in R^4
% from the flat outputs and their time derivatives, without integration.
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
%   x(4:6)   = v_b         Body velocities via R^T * v_I.
%   x(7:8)   = phi, theta  Derived from thrust vector (z_B) and desired heading.
%   x(9)     = psi         Euler yaw extracted from R (psi_check).
%   x(10:12) = p, q, r     Body rates derived from position jerk and yaw rate.
%
% Input reconstruction u = i(sigma, ..., ddddsigma)
%   u(1) = T   from ||F_net||
%   u(2) = L   from snap + Euler eq.
%   u(3) = M   from snap + Euler eq.
%   u(4) = N   from snap + Euler eq.
%
% Singularities
%   T = 0 (free fall): thrust direction undefined.
%   z_B parallel to x_C: gimbal-like singularity at 90° pitch aligned with yaw.

    methods (Static)

            function [xref, uref] = map(mp, sigma, dsigma, ddsigma, dddsigma, ddddsigma)
            % MAP  Algebraic reconstruction of reference state and input.
            %
            %   [xref, uref] = FlatnessMap.map(mp, sigma, dsigma, ddsigma, dddsigma, ddddsigma)
            %
            %   Inputs
            %     mp         : ModelParameters struct (fields m,g,Jxx,Jyy,Jzz)
            %     sigma      : [4x1] flat outputs     [x_I; y_I; z_I; psi]
            %     dsigma     : [4x1] 1st derivatives  [vx_I; vy_I; vz_I; psi_dot] 
            %     ddsigma    : [4x1] 2nd derivatives  [ax_I; ay_I; az_I; psi_ddot]
            %     dddsigma   : [4x1] 3rd derivatives  [jx_I; jy_I; jz_I; psi_dddot] (jerk)
            %     ddddsigma  : [4x1] 4th derivatives  [sx_I; sy_I; sz_I; psi_ddddot] (snap)
            %
            %   Outputs
            %     xref : [12x1] reference state  [p_I; v_b; phi; theta; psi; p; q; r]
            %     uref : [4x1]  reference input  [T; L; M; N]

            %% Parameters
            g  = mp.g;
            m  = mp.m;

            Jxx = mp.Jxx; 
            Jyy = mp.Jyy; 
            Jzz = mp.Jzz;
            J = diag([Jxx Jyy Jzz]);


            %% Flat Outputs
            psi     = sigma(4); % psi [rad]
            psi_dot = dsigma(4); % dpsi [rad/s]
            psi_dd  = ddsigma(4); % ddpsi [rad/s^2]
                
            p_I = sigma(1:3);     % inertial position [m]
            v_I = dsigma(1:3);    % inertial velocity [m/s]
            a_I = ddsigma(1:3);   % inertial acceleration [m/s^2]
            j_I = dddsigma(1:3);  % inertial jerk [m/s^3]
            s_I = ddddsigma(1:3); % inertial snap [m/s^4]


             %%  Net force vector in NED frame
             %
             %  From Newton in inertial frame:
             %
             %  m * a_I = m * g *e_z + R * [0;0;-T]
             %
             %  Solving for the thrust vector:
             %
             %    -T * z_B^I  =  m*(a_I - [0;0;g]) 
             %                =: F_net
             %
             % F_net = m * a_I - m * g * e_z 
             %
             %  Note: g appears with + sign because NED convention:
             %        e_z = [0;0;1] points Down, and gravity IS m*g*e_z.

             F_net = m * (a_I - [0; 0; g]);
             T_ref = norm(F_net);

             % Check for free fall
             assert(T_ref > 1e-4, ...
                 'FlatnessMap:ZeroThrust', ...
                 'T = %.4f N: trajectory passes through free fall (T=0). Increase thrust margin.', T_ref);

             z_B = -F_net / T_ref;   % unit vector: z_B axis in NED frame

                
             %% Euler Angles [phi, theta, psi]
             %
             %  psi is the 4th flat output (known exactly).
             %
             %  phi and theta are derived from z_B^I and the desired yaw,
             %
             %  1. x_C = [cos(psi), sin(psi), 0]  (horizontal "heading" vector)
             %  2. y_B^I = (z_B^I × x_C) / ||(z_B^I × x_C) || 
             %  3. x_B^I = y_B^I × z_B^I
             %  4. R = [x_B^I | y_B^I | z_B^I]   (columns)
             %  5. Extract phi, theta, psi_check from R (Z-Y-X decomposition)
             %
             %  The singularity occurs when z_B^I || x_C, i.e. when the
             %  drone pitches ±90° aligned with the yaw direction.
                
             % Horizontal heading vector
             x_C = [cos(psi); sin(psi); 0];
             y_B_cross = cross(z_B, x_C);
             y_B_norm  = norm(y_B_cross);

             % Check for Attitude Singularity
             assert(y_B_norm > 1e-6, ...
                 'FlatnessMap:AttitudeSingularity', ...
                 'z_B parallel to x_C: pitch singularity near 90 deg aligned with yaw.');

              y_B = y_B_cross / y_B_norm;
              x_B = cross(y_B, z_B);

              % Rotation matrix  R = [x_B | y_B | z_B]  (body axes in NED)
              R = [x_B, y_B, z_B];

              % Extract Euler angles  (Z-Y-X convention)
              %   R(3,1) = -sin(theta)
              %   R(3,2) =  sin(phi)*cos(theta)
              %   R(3,3) =  cos(phi)*cos(theta)
              psi_check = atan2(R(2,1), R(1,1));
              theta = asin(-R(3,1));
              phi   = atan2(R(3,2), R(3,3));
                
              % Yaw error between psi (flat ouptut) and reconstructed psi
              % (psi_check)
              yaw_err = atan2( ...
                  sin(psi - psi_check), ...
                  cos(psi - psi_check));
                
              % Reconstruction Error
              %assert(abs(yaw_err) < 1e-3, ...
                  %'FlatnessMap:YawMismatch', ...
                  %'Yaw reconstruction inconsistency.');

              cthe = cos(theta);

              % Check for Gimbal Lock
              assert(abs(cthe) > 1e-6, ...
                 'FlatnessMap:GimbalLock', ...
                 'cos(theta) ~ 0.');
              

              %% Body frame velocities  v_b = R^T * v_I
              %
              %  The state stores velocities in body FRD frame.
              %  Given v_I (from flat output derivatives) and R:
              %
              %    v_b = R^T * v_I = [x_B^T; y_B^T; z_B^T] * v_I
              %
              %  This is the exact inversion of:  v_I = R * v_b  
              v_b = R' * v_I;   % [xb_dot; yb_dot; zb_dot]


              %% Angular accelerations  [p_dot, q_dot, r_dot]
              %
              %  The angular acceleration is reconstructed by differentiating the
              %  rotation kinematics:
              %
              %      omega_hat = R^T * R_dot
              %
              %  Differentiating:
              %
              %      omega_hat_dot = R^T * R_ddot - omega_hat^2
              %
              %  To compute R_ddot, differentiate the body axes again.
              %
              %  Since:
              %
              %      F_ddot = m * s_I
              %
              %  and:
              %
              %      T_ddot = -(z_B_dot^T * F_dot + z_B^T * F_ddot)
              %
              %  then:
              %
              %      z_B_ddot = -(F_ddot + T_ddot*z_B + 2*T_dot*z_B_dot)/T
              %
              %  The commanded heading acceleration is:
              %
              %      x_C_ddot =
              %          psi_dd * [-sin(psi); cos(psi); 0]
              %          + psi_dot^2 * [-cos(psi); -sin(psi); 0]
              %
              %  Define:
              %
              %      c_ddot =
              %          z_B_ddot x x_C
              %          + 2*(z_B_dot x x_C_dot)
              %          + z_B x x_C_ddot
              %
              %  The second derivative of the normalized vector y_B is:
              %
              %      y_B_ddot = P_y * c_ddot / norm(c_vec) - ( ...
              %                 2 * y_B_dot * (y_B' * c_dot) ...
              %                 + y_B * (y_B_dot' * c_dot) ...
              %                 ) / norm(c_vec);
              %
              %  Then:
              %
              %      x_B_ddot =
              %          y_B_ddot x z_B
              %          + 2*(y_B_dot x z_B_dot)
              %          + y_B x z_B_ddot
              %
              %  Therefore:
              %
              %      R_ddot = [x_B_ddot, y_B_ddot, z_B_ddot]
              %
              %  and:
              %
              %      omega_hat_dot = R^T * R_ddot - omega_hat^2
              %
              %  Extract angular accelerations:
              %
              %      p_dot = omega_hat_dot(3,2)
              %      q_dot = omega_hat_dot(1,3)
              %      r_dot = omega_hat_dot(2,1)

              F_dot = m* j_I;
              T_dot = -(z_B' * F_dot);

              z_B_dot = -(F_dot + T_dot*z_B) / T_ref;
              x_C_dot = psi_dot * [-sin(psi); cos(psi); 0];

              c_vec = cross(z_B, x_C);
              c_dot = cross(z_B_dot, x_C) + cross(z_B, x_C_dot);

              y_B_dot = (eye(3) - y_B*y_B') * c_dot / norm(c_vec);

              % x_B = y_B x z_B
              x_B_dot = cross(y_B_dot, z_B) + cross(y_B, z_B_dot);

              % R DOT
              R_dot = [x_B_dot y_B_dot z_B_dot];

              % omega_hat = R' * R_dot
              omega_hat = R' * R_dot;

              % extract body rates
              p_ref = omega_hat(3,2);
              q_ref = omega_hat(1,3);
              r_ref = omega_hat(2,1);
                
              omega = [p_ref; q_ref; r_ref];

              F_ddot = m * s_I;
              T_ddot = -(z_B_dot' * F_dot + z_B' * F_ddot);
                
              z_B_ddot = -(F_ddot + T_ddot*z_B + 2*T_dot*z_B_dot) / T_ref;
              x_C_ddot = psi_dd * [-sin(psi); cos(psi); 0] + psi_dot^2 * [-cos(psi); -sin(psi); 0];
              c_ddot = cross(z_B_ddot, x_C) + 2*cross(z_B_dot, x_C_dot) + cross(z_B, x_C_ddot);

              % normalized vector second derivative
              P_y = eye(3) - y_B*y_B';

              % complete:
              y_B_ddot = P_y * c_ddot / norm(c_vec) - ( ...
                  2 * y_B_dot * (y_B' * c_dot) ...
                  + y_B * (y_B_dot' * c_dot) ...
                  ) / norm(c_vec);

              x_B_ddot = ...
                  cross(y_B_ddot, z_B) ...
                  + 2*cross(y_B_dot, z_B_dot) ...
                  + cross(y_B, z_B_ddot);

              % R DDOT
              R_ddot = [x_B_ddot y_B_ddot z_B_ddot];

              % Differentiate:
              %   omega_hat = R'R_dot
              %   omega_hat_dot = R'R_ddot - omega_hat^2
              omega_hat_dot = R' * R_ddot - omega_hat * omega_hat;

              p_dot_ref = omega_hat_dot(3,2);
              q_dot_ref = omega_hat_dot(1,3);
              r_dot_ref = omega_hat_dot(2,1);

              omega_dot = [p_dot_ref; q_dot_ref; r_dot_ref];


              %% Torques  [L, M, N]
              %
              %  Control torques are reconstructed from the rigid-body rotational
              %  dynamics:
              %
              %      J*omega_dot + omega x (J*omega) = tau
              %
              %  where:
              %
              %      omega     = [p; q; r]
              %      omega_dot = [p_dot; q_dot; r_dot]
              %
              %  and:
              %
              %      J = diag(Jxx, Jyy, Jzz)
              %
              %  Therefore:
              %
              %      tau = J*omega_dot + omega x (J*omega)
              %
              %  whose components are:
              %
              %      L = tau(1)
              %      M = tau(2)
              %      N = tau(3)
              %
              %  This reconstruction is exact and directly obtained from the desired
              %  snap trajectory and attitude evolution.

              tau = J*omega_dot + cross(omega, J*omega);

              L_ref = tau(1);
              M_ref = tau(2);
              N_ref = tau(3);


              %% ASSEMBLE STATE AND INPUT REFERENCE VECTORS
        
              %  State vector (12x1) 
              %  [X_I, Y_I, Z_I, xb_dot, yb_dot, zb_dot, phi, theta, psi, p, q, r]
              xref = [p_I;              % 1-3  position NED
                      v_b;              % 4-6  body velocities
                      phi; theta; psi_check;  % 7-9  Euler angles
                      p_ref; q_ref; r_ref];  % 10-12 body rates

              %  Input vector (4x1):
              %  [T, L, M, N]
              uref = [T_ref; L_ref; M_ref; N_ref];

            end
    end
end