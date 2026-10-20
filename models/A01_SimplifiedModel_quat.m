classdef A01_SimplifiedModel_quat
    % 13-state rigid body quadrotor — quaternion attitude representation.
    % Quaternion replaces Z-Y-X Euler angles, NO gimbal lock.
    %
    % Reference frames
    %   Inertial : NED  (North-East-Down)
    %   Body     : FRD  (Forward-Right-Down)
    %   Rotation : q = [q0; q1; q2; q3], scalar first
    %              q represents the rotation body -> NED: v_I = R(q) * v_b
    %
    % State vector  x in R^13
    %   x(1)  = X_I    [m]      North position  (NED)
    %   x(2)  = Y_I    [m]      East position  (NED)
    %   x(3)  = Z_I    [m]      Down  position  (NED, positive downward)
    %   x(4)  = xb_dot [m/s]    forward velocity (FRD)
    %   x(5)  = yb_bot [m/s]    rightward velocity (FRD)
    %   x(6)  = zb_dot [m/s]    downward velocity (FRD)
    %   x(7)  = q0     [-]      quaternion scalar
    %   x(8)  = q1     [-]      quaternion vector (x)
    %   x(9)  = q2     [-]      quaternion vector (y)
    %   x(10) = q3     [-]      quaternion vector (z)
    %   x(11) = p      [rad/s]  roll rate (body)
    %   x(12) = q_r    [rad/s]  pitch rate (body)
    %   x(13) = r      [rad/s]  yaw rate (body)
    %
    % Input vector  u in R^4  
    %   u(1) = T  [N]    collective thrust along -z_B
    %   u(2) = L  [Nm]   roll torque
    %   u(3) = M  [Nm]   pitch torque
    %   u(4) = N  [Nm]   yaw torque
    %
    % Hover trim:
    %   x* = [0;0;0; 0;0;0; 1;0;0;0; 0;0;0],  u* = [m*g; 0; 0; 0]

    properties
        mp
    end

    methods

        function obj = A01_SimplifiedModel_quat(params)
            obj.mp = params;
        end

        function dxdt = dynamics(obj, t, x, u)
            % DYNAMICS  Evaluate dxdt = f(x, u).
            %
            %   Inputs
            %     t  : scalar time [s]
            %     x  : 13x1 state vector
            %     u  : 4x1  input [T; L; M; N]
            %   Output
            %     dxdt : 13x1

            % state
            xb_dot = x(4);  yb_dot = x(5);  zb_dot = x(6);
            q0  = x(7);  q1  = x(8);  q2  = x(9);  q3  = x(10);
            p   = x(11); q_r = x(12); r   = x(13);   % q_r = pitch rate

            % input
            T = u(1);  L = u(2);  M = u(3);  N = u(4);

            % parameters
            m   = obj.mp.m;   g   = obj.mp.g;
            Jxx = obj.mp.Jxx; Jyy = obj.mp.Jyy; Jzz = obj.mp.Jzz;


            %% BLOCK 1 — Position kinematics:  p_I_dot = R(q) * v_b
            %
            % Rotation matrix R (body FRD -> NED), from quaternion:
            %   R = (q0^2 - ||q_v||^2)*I + 2*q_v*q_v' + 2*q0*S(q_v)
            %
            % Rows of R applied to [xb_dot; yb_dot; zb_dot]:

            dX_I = (q0^2+q1^2-q2^2-q3^2)*xb_dot + 2*(q1*q2-q0*q3)*yb_dot + 2*(q1*q3+q0*q2)*zb_dot;
            dY_I = 2*(q1*q2+q0*q3)*xb_dot + (q0^2-q1^2+q2^2-q3^2)*yb_dot + 2*(q2*q3-q0*q1)*zb_dot;
            dZ_I = 2*(q1*q3-q0*q2)*xb_dot + 2*(q2*q3+q0*q1)*yb_dot + (q0^2-q1^2-q2^2+q3^2)*zb_dot;


            %% BLOCK 2 — Translational dynamics (Newton in body frame)
            %
            % Gravity in body frame: g_b = R(q)^T * [0; 0; g]
            %   g_b(1) = 2*(q1*q3 - q0*q2)*g   [= -g*sin(theta)]
            %   g_b(2) = 2*(q2*q3 + q0*q1)*g   [= +g*sin(phi)*cos(theta)]
            %   g_b(3) = (q0^2-q1^2-q2^2+q3^2)*g  [= +g*cos(phi)*cos(theta)]
            %
            % omega x v_b:
            %   xb_ddot = r*yb_dot - q_r*zb_dot + g_b(1)
            %   yb_ddot = p*zb_dot - r*xb_dot  + g_b(2)
            %   zb_ddot = q_r*xb_dot - p*yb_dot + g_b(3) - T/m

            g_b1 = 2*(q1*q3 - q0*q2)*g;
            g_b2 = 2*(q2*q3 + q0*q1)*g;
            g_b3 = (q0^2 - q1^2 - q2^2 + q3^2)*g;

            dxb_dot = r*yb_dot   - q_r*zb_dot + g_b1;
            dyb_dot = p*zb_dot   - r*xb_dot   + g_b2;
            dzb_dot = q_r*xb_dot - p*yb_dot   + g_b3 - T/m;


            %% BLOCK 3 — Quaternion kinematics
            %
            % q_dot = (1/2) * q x [0; omega_b]  (Hamilton product)
            %
            %   dq0 = -(1/2)*(q1*p + q2*q_r + q3*r)
            %   dq1 =  (1/2)*(q0*p + q2*r  - q3*q_r)
            %   dq2 =  (1/2)*(q0*q_r + q3*p - q1*r)
            %   dq3 =  (1/2)*(q0*r  + q1*q_r - q2*p)
            %
            % No gimbal lock; 

            dq0 = 0.5 * (-(q1*p + q2*q_r + q3*r));
            dq1 = 0.5 * ( q0*p + q2*r   - q3*q_r);
            dq2 = 0.5 * ( q0*q_r + q3*p - q1*r  );
            dq3 = 0.5 * ( q0*r  + q1*q_r - q2*p );


            %% BLOCK 4 — Rotational dynamics (Euler equations)  

            dp  = (1/Jxx) * (L   - (Jzz-Jyy)*q_r*r);
            dqr = (1/Jyy) * (M - (Jxx-Jzz)*p*r  );
            dr  = (1/Jzz) * (N - (Jyy-Jxx)*p*q_r);


            %% Assemble
            dxdt = [dX_I; dY_I; dZ_I;
                    dxb_dot; dyb_dot; dzb_dot;
                    dq0;  dq1;  dq2;  dq3;
                    dp;   dqr;  dr];
        end

        function [A, B] = get_jacobians(obj, x_lin, u_lin)
            % GET_JACOBIANS  Analytical linearization df/dx (13x13) and df/du (13x4).
            %
            % At hover (v_b=0, q=[1;0;0;0], omega=0, T=m*g):
            %   A(1,4)=1  A(2,5)=1  A(3,6)=1        (kinematic: R=I)
            %   A(4,9)=-2g  A(5,7)=2g  A(6,7)=2g    (quaternion gravity coupling)
            %   B(6,1)=-1/m  B(11,2)=1/Jxx  B(12,3)=1/Jyy  B(13,4)=1/Jzz

            % unpack
            xb_dot = x_lin(4);  yb_dot = x_lin(5);  zb_dot = x_lin(6);
            q0  = x_lin(7);  q1  = x_lin(8);  q2  = x_lin(9);  q3  = x_lin(10);
            p   = x_lin(11); q_r = x_lin(12); r   = x_lin(13);

            % parameters
            m   = obj.mp.m;   g   = obj.mp.g;
            Jxx = obj.mp.Jxx; Jyy = obj.mp.Jyy; Jzz = obj.mp.Jzz;
            
            % initialize matrices
            A = zeros(13, 13);
            B = zeros(13, 4);


            %% Block 1 rows (1-3): P_I_dot = R(q) * v_b
            % Part (a): d(R*vb)/d(vb) = R  [cols 4,5,6]
            A(1,4) = q0^2+q1^2-q2^2-q3^2;
            A(1,5) = 2*(q1*q2-q0*q3);
            A(1,6) = 2*(q1*q3+q0*q2);

            A(2,4) = 2*(q1*q2+q0*q3);
            A(2,5) = q0^2-q1^2+q2^2-q3^2;
            A(2,6) = 2*(q2*q3-q0*q1);

            A(3,4) = 2*(q1*q3-q0*q2);
            A(3,5) = 2*(q2*q3+q0*q1);
            A(3,6) = q0^2-q1^2-q2^2+q3^2;

            % Part (b): d(R*v_b)/d(q)  [cols 7,8,9,10]
            A(1,7) = 2*( q0*xb_dot - q3*yb_dot + q2*zb_dot);
            A(1,8) = 2*( q1*xb_dot + q2*yb_dot + q3*zb_dot);
            A(1,9) = 2*(-q2*xb_dot + q1*yb_dot + q0*zb_dot);
            A(1,10)= 2*(-q3*xb_dot - q0*yb_dot + q1*zb_dot);

            A(2,7) = 2*( q3*xb_dot + q0*yb_dot - q1*zb_dot);
            A(2,8) = 2*( q2*xb_dot - q1*yb_dot - q0*zb_dot);
            A(2,9) = 2*( q1*xb_dot + q2*yb_dot + q3*zb_dot);
            A(2,10)= 2*( q0*xb_dot - q3*yb_dot + q2*zb_dot);

            A(3,7) = 2*(-q2*xb_dot + q1*yb_dot + q0*zb_dot);
            A(3,8) = 2*( q3*xb_dot + q0*yb_dot - q1*zb_dot);
            A(3,9) = 2*(-q0*xb_dot + q3*yb_dot - q2*zb_dot);
            A(3,10)= 2*( q1*xb_dot + q2*yb_dot + q3*zb_dot);


            %% Block 2 rows (4-6): translational dynamics
            % Coriolis velocity cross-coupling
            A(4,5) =  r;     A(4,6) = -q_r;
            A(5,4) = -r;     A(5,6) =  p;
            A(6,4) =  q_r;   A(6,5) = -p;

            % Gravity-quaternion coupling: d(g_b)/d(q)  [cols 7,8,9,10]
            %   g_b1 = 2*(q1*q3-q0*q2)*g
            A(4,7) = -2*q2*g;  A(4,8) =  2*q3*g;  A(4,9) = -2*q0*g;  A(4,10)=  2*q1*g;
            %   g_b2 = 2*(q2*q3+q0*q1)*g
            A(5,7) =  2*q1*g;  A(5,8) =  2*q0*g;  A(5,9) =  2*q3*g;  A(5,10)=  2*q2*g;
            %   g_b3 = (q0^2-q1^2-q2^2+q3^2)*g
            A(6,7) =  2*q0*g;  A(6,8) = -2*q1*g;  A(6,9) = -2*q2*g;  A(6,10)=  2*q3*g;

            % Coriolis rate coupling  [cols 11,12,13]
            A(4,12) = -zb_dot;   A(4,13) =  yb_dot;
            A(5,11) =  zb_dot;   A(5,13) = -xb_dot;
            A(6,11) = -yb_dot;   A(6,12) =  xb_dot;


            %% Block 3 rows (7-10): quaternion kinematics
            % d(q_dot)/d(q) = (1/2)*Omega(omega)  [cols 7,8,9,10]
            A(7,7)  = 0;      A(7,8)  = -p/2;   A(7,9)  = -q_r/2; A(7,10) = -r/2;
            A(8,7)  =  p/2;   A(8,8)  = 0;      A(8,9)  =  r/2;   A(8,10) = -q_r/2;
            A(9,7)  =  q_r/2; A(9,8)  = -r/2;   A(9,9)  = 0;      A(9,10) =  p/2;
            A(10,7) =  r/2;   A(10,8) =  q_r/2; A(10,9) = -p/2;   A(10,10)= 0;

            % d(q_dot)/d(omega) = (1/2)*Xi(q)  [cols 11,12,13]
            A(7,11) = -q1/2;  A(7,12)  = -q2/2;  A(7,13)  = -q3/2;
            A(8,11) =  q0/2;  A(8,12)  = -q3/2;  A(8,13)  =  q2/2;
            A(9,11) =  q3/2;  A(9,12)  =  q0/2;  A(9,13)  = -q1/2;
            A(10,11)= -q2/2;  A(10,12) =  q1/2;  A(10,13) =  q0/2;


            %% Block 4 rows (11-13): rotational dynamics (identical to A00, shifted)
            A(11,12) = -(Jzz-Jyy)*r/Jxx;    A(11,13) = -(Jzz-Jyy)*q_r/Jxx;
            A(12,11) = -(Jxx-Jzz)*r/Jyy;    A(12,13) = -(Jxx-Jzz)*p/Jyy;
            A(13,11) = -(Jyy-Jxx)*q_r/Jzz;  A(13,12) = -(Jyy-Jxx)*p/Jzz;


            %% B matrix: df/du
            B(6,1)  = -1/m;
            B(11,2) =  1/Jxx;
            B(12,3) =  1/Jyy;
            B(13,4) =  1/Jzz;
        end

    end
end