classdef A02_SimplifiedModel_Deployment
    % SIMPLIFIED MODEL - only rigid body  Nonlinear 6-DoF quadrotor.
    %
    % Newton-Euler physics, velocities are expressed in the inertial NED frame.
    % Same physics as A00_SimplifiedModel: only the velocity state changes
    % (DroneCmd getVelocity: NED).
    %
    % Reference frames
    %   Inertial : NED  (North-East-Down), fixed to ground
    %   Body     : FRD  (Forward-Right-Down), attached to CoG
    %
    % State vector  x in R^12  (column)
    %   x(1)  = X_I    [m]      North position  (NED)
    %   x(2)  = Y_I    [m]      East  position  (NED)
    %   x(3)  = Z_I    [m]      Down  position  (NED, positive downward)
    %   x(4)  = v_x    [m/s]    North velocity  (inertial NED)
    %   x(5)  = v_y    [m/s]    East  velocity  (inertial NED)
    %   x(6)  = v_z    [m/s]    Down  velocity  (inertial NED)
    %   x(7)  = phi    [rad]    Roll  angle (Z-Y-X Euler)
    %   x(8)  = theta  [rad]    Pitch angle (positive nose-up)
    %   x(9)  = psi    [rad]    Yaw   angle
    %   x(10) = p      [rad/s]  Roll  rate (body)
    %   x(11) = q      [rad/s]  Pitch rate (body)
    %   x(12) = r      [rad/s]  Yaw   rate (body)
    %
    % Input vector  u in R^4  (column)
    %   u(1)  = T   [N]     Collective thrust along -z_B
    %   u(2)  = L   [N*m]   Roll  torque about x_B
    %   u(3)  = M   [N*m]   Pitch torque about y_B
    %   u(4)  = N   [N*m]   Yaw   torque about z_B
    %
    % Hover:
    %   x* = zeros(12,1),   u* = [m*g; 0; 0; 0]

    properties
        mp
    end

    methods

        function obj = A02_SimplifiedModel_Deployment(params)
            % Constructor.
            %   params : ModelParameters object
            obj.mp = params;
        end

        function dxdt = dynamics(obj, t, x, u)
            % DYNAMICS  Evaluate the nonlinear ODE  dxdt = f(x,u).
            %
            %   Inputs
            %     t  : scalar time [s]
            %     x  : 12x1 state vector
            %     u  : 4x1  input [T; L; M; N]
            %   Output
            %     dxdt : 12x1

            % state
            v_x   = x(4);   % North velocity (NED)
            v_y   = x(5);   % East  velocity (NED)
            v_z   = x(6);   % Down  velocity (NED)
            phi   = x(7);
            theta = x(8);
            psi   = x(9);
            p     = x(10);
            q     = x(11);
            r     = x(12);

            % input
            T = u(1);
            L = u(2);
            M = u(3);
            N = u(4);

            % parameters
            m   = obj.mp.m;
            g   = obj.mp.g;
            Jxx = obj.mp.Jxx;
            Jyy = obj.mp.Jyy;
            Jzz = obj.mp.Jzz;

            % trigonometric functions
            sphi = sin(phi);   cphi = cos(phi);
            sthe = sin(theta); cthe = cos(theta);
            spsi = sin(psi);   cpsi = cos(psi);

            % gimbal lock assertion
            if isnumeric(theta)
                assert(abs(cthe) > 1e-6, ...
                    'QuadrotorModel:GimbalLock', ...
                    'theta = %.2f rad: singular (gimbal lock). Use 13-state model.', theta);
            end

            %% BLOCK 1 — Position kinematics
            %
            %   P_I_dot = v_I
            %
            % Velocities are already expressed in NED, trivial

            dX_I = v_x;
            dY_I = v_y;
            dZ_I = v_z;


            %% BLOCK 2 — Translational dynamics (Newton in inertial frame)
            %
            %   m*(dv_I/dt) = m*g*e_z + R(phi,theta,psi) * [0; 0; -T]
            %
            %   R(:,3) = [cphi*sthe*cpsi + sphi*spsi]
            %            [cphi*sthe*spsi - sphi*cpsi]
            %            [cphi*cthe                 ]
            %
            %   dv_I = g*e_z - (T/m)*R(:,3)
            %     v_x_d = -(T/m)*(cphi*sthe*cpsi + sphi*spsi)
            %     v_y_d = -(T/m)*(cphi*sthe*spsi - sphi*cpsi)
            %     v_z_d =  g - (T/m)*cphi*cthe

            dv_x = -(T/m) * (cphi*sthe*cpsi + sphi*spsi);
            dv_y = -(T/m) * (cphi*sthe*spsi - sphi*cpsi);
            dv_z =  g - (T/m) * cphi*cthe;


            %% BLOCK 3 — Rotational kinematics (Euler Z-Y-X)

            tthe = sthe / cthe;
            iche = 1.0  / cthe;

            dphi   = p + sphi*tthe*q + cphi*tthe*r;
            dtheta =     cphi*q      - sphi*r;
            dpsi   =     sphi*iche*q + cphi*iche*r;


            %% BLOCK 4 — Rotational dynamics (Euler equations)

            dp_ = (1/Jxx) * (L - (Jzz - Jyy)*q*r);
            dq_ = (1/Jyy) * (M - (Jxx - Jzz)*p*r);
            dr_ = (1/Jzz) * (N - (Jyy - Jxx)*p*q);


            %% assemble
            dxdt = [dX_I; dY_I; dZ_I;
                    dv_x; dv_y; dv_z;
                    dphi; dtheta; dpsi;
                    dp_; dq_; dr_];
        end


        function [A, B] = get_jacobians(obj, x_lin, u_lin)
            % GET_JACOBIANS  Analytical linearization at (x_lin, u_lin).
            %
            %   Returns A = df/dx (12x12) and B = df/du (12x4).
            %
            %
            %   At hover (phi=theta=psi=0, p=q=r=0, T=mg, v_I=0):
            %     A(1,4)  =  1         (X_I_dot / v_x)
            %     A(2,5)  =  1         (Y_I_dot / v_y)
            %     A(3,6)  =  1         (Z_I_dot / v_z)
            %     A(4,8)  = -g         (v_x_dot / theta: thrust tilt)
            %     A(5,7)  = +g         (v_y_dot / phi:   thrust tilt)
            %     B(6,1)  = -1/m
            %     B(10,2) = 1/Jxx   B(11,3)=1/Jyy   B(12,4)=1/Jzz

            % unpack
            phi   = x_lin(7);
            theta = x_lin(8);
            psi   = x_lin(9);
            p     = x_lin(10);
            q     = x_lin(11);
            r     = x_lin(12);
            T     = u_lin(1);

            % parameters
            m   = obj.mp.m;
            Jxx = obj.mp.Jxx;
            Jyy = obj.mp.Jyy;
            Jzz = obj.mp.Jzz;

            % trigonometric functions
            sphi = sin(phi);   cphi = cos(phi);
            sthe = sin(theta); cthe = cos(theta);
            spsi = sin(psi);   cpsi = cos(psi);
            tthe = sthe / cthe;
            iche = 1.0  / cthe;

            % gimbal lock
            assert(abs(cthe) > 1e-6, ...
                'QuadrotorModel:GimbalLock', ...
                'theta = %.2f rad: singular (gimbal lock). Use 13-state model.', theta);

            % initialize matrices
            A = zeros(12, 12);
            B = zeros(12, 4);

            %% Matrix A
            %% BLOCK 1 rows (1-3): P_I_dot = v_I

            A(1,4) = 1;
            A(2,5) = 1;
            A(3,6) = 1;


            %% BLOCK 2 rows (4-6): translational dynamics in inertial frame
            %
            %   dv_x = -(T/m)*(cphi*sthe*cpsi + sphi*spsi)
            %   dv_y = -(T/m)*(cphi*sthe*spsi - sphi*cpsi)
            %   dv_z =  g - (T/m)*cphi*cthe
            %
            % Derivatives w.r.t. angles (thrust direction):
            %   d(dv_x)/d(phi)   =  (T/m)*(cpsi*sthe*sphi - spsi*cphi)
            %   d(dv_x)/d(theta) = -(T/m)*cpsi*cthe*cphi
            %   d(dv_x)/d(psi)   =  (T/m)*(spsi*sthe*cphi - cpsi*sphi)
            %   d(dv_y)/d(phi)   =  (T/m)*(spsi*sthe*sphi + cpsi*cphi)
            %   d(dv_y)/d(theta) = -(T/m)*spsi*cthe*cphi
            %   d(dv_y)/d(psi)   = -(T/m)*(cpsi*sthe*cphi + spsi*sphi)
            %   d(dv_z)/d(phi)   =  (T/m)*cthe*sphi
            %   d(dv_z)/d(theta) =  (T/m)*sthe*cphi
            %   d(dv_z)/d(psi)   =  0

            % row 4: dv_x
            A(4,7) =  (T/m) * (cpsi*sthe*sphi - spsi*cphi);
            A(4,8) = -(T/m) * cpsi*cthe*cphi;
            A(4,9) =  (T/m) * (spsi*sthe*cphi - cpsi*sphi);

            % row 5: dv_y
            A(5,7) =  (T/m) * (spsi*sthe*sphi + cpsi*cphi);
            A(5,8) = -(T/m) * spsi*cthe*cphi;
            A(5,9) = -(T/m) * (cpsi*sthe*cphi + spsi*sphi);

            % row 6: dv_z
            A(6,7) =  (T/m) * cthe*sphi;
            A(6,8) =  (T/m) * sthe*cphi;
            % A(6,9) = 0; 


            %% A — BLOCK 3 rows (7-9): Euler kinematics

            A(7,7)  = cphi*tthe*q - sphi*tthe*r;
            A(7,8)  = sphi*(iche^2)*q + cphi*(iche^2)*r;
            A(7,10) = 1;
            A(7,11) = sphi*tthe;
            A(7,12) = cphi*tthe;

            A(8,7)  = -sphi*q - cphi*r;
            A(8,11) = cphi;
            A(8,12) = -sphi;

            A(9,7)  = cphi*iche*q - sphi*iche*r;
            A(9,8)  = sphi*sthe*(iche^2)*q + cphi*sthe*(iche^2)*r;
            A(9,11) = sphi*iche;
            A(9,12) = cphi*iche;


            %% BLOCK 4 rows (10-12): rotational dynamics

            % A(10,10) = 0;
            A(10,11) = -(Jzz-Jyy)*r/Jxx;
            A(10,12) = -(Jzz-Jyy)*q/Jxx;

            % A(11,11) = 0;
            A(11,10) = -(Jxx-Jzz)*r/Jyy;
            A(11,12) = -(Jxx-Jzz)*p/Jyy;

            % A(12,12) = 0;
            A(12,10) = -(Jyy-Jxx)*q/Jzz;
            A(12,11) = -(Jyy-Jxx)*p/Jzz;


            %% Matrix B
            %% Input Jacobian df/du
            %
            %  Thrust T acts on all three inertial accelerations through R(:,3):
            %    d(dv_I)/dT = -(1/m) * R(:,3)

            B(4,1)  = -(1/m) * (cphi*sthe*cpsi + sphi*spsi);
            B(5,1)  = -(1/m) * (cphi*sthe*spsi - sphi*cpsi);
            B(6,1)  = -(1/m) * cphi*cthe;
            B(10,2) = 1/Jxx;
            B(11,3) = 1/Jyy;
            B(12,4) = 1/Jzz;

        end
    end
end