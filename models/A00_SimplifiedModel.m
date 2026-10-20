classdef A00_SimplifiedModel
    % SIMPLIFIED MODEL - only rigid body  Nonlinear 6-DoF quadrotor (body velocities model).
    %
    % Newton-Euler physics, velocities
    % are expressed in the body FRD frame.
    %
    % Reference frames
    %   Inertial : NED  (North-East-Down), fixed to ground
    %   Body     : FRD  (Forward-Right-Down), attached to CoG
    %
    % State vector  x in R^12  (column)
    %   x(1)  = X_I    [m]      North position  (NED)
    %   x(2)  = Y_I    [m]      East  position  (NED)
    %   x(3)  = Z_I    [m]      Down  position  (NED, positive downward)
    %   x(4)  = xb_dot [m/s]    Forward  velocity (body FRD)
    %   x(5)  = yb_dot [m/s]    Rightward velocity (body FRD)
    %   x(6)  = zb_dot [m/s]    Downward velocity  (body FRD)
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
    % Hover trim:
    %   x* = zeros(12,1),   u* = [m*g; 0; 0; 0]

    properties
        mp
    end

    methods

        function obj = A00_SimplifiedModel(params)
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
            xb_d  = x(4);   % xb_dot: body x velocity
            yb_d  = x(5);   % yb_dot: body y velocity
            zb_d  = x(6);   % zb_dot: body z velocity
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

            % gimbal lock 
            assert(abs(cthe) > 1e-6, ...
                'QuadrotorModel:GimbalLock', ...
                'theta = %.2f rad: singular (gimbal lock). Use 13-state model.', theta);

            %% BLOCK 1 — Position kinematics
            %
            %   P_I_dot = R(phi,theta,psi) * v_b
            %
            % Full matrix (body FRD -> NED), rotation sequence Z-Y-X:
            %
            %   R = [cthe*cpsi,  sphi*sthe*cpsi - cphi*spsi,  cphi*sthe*cpsi + sphi*spsi]
            %       [cthe*spsi,  sphi*sthe*spsi + cphi*cpsi,  cphi*sthe*spsi - sphi*cpsi]
            %       [-sthe,      sphi*cthe,                   cphi*cthe                 ]
            %
            % Each row of R applied to [xb_d; yb_d; zb_d]:

            dX_I = cthe*cpsi*xb_d + (sphi*sthe*cpsi - cphi*spsi)*yb_d + (cphi*sthe*cpsi + sphi*spsi)*zb_d;
            dY_I = cthe*spsi*xb_d + (sphi*sthe*spsi + cphi*cpsi)*yb_d + (cphi*sthe*spsi - sphi*cpsi)*zb_d;
            dZ_I = -sthe*xb_d     + sphi*cthe*yb_d                    + cphi*cthe*zb_d;


            %% BLOCK 2 — Translational dynamics (Newton in body frame)
            %
            %   m*(dv_b/dt) = F_b - m*(omega_b x v_b)
            %
            %   F_b = gravity_in_body + thrust 
            %       = [ -g*sin(theta)         ] + [   0   ] 
            %         [  g*sin(phi)*cos(theta) ]   [   0   ]   
            %         [  g*cos(phi)*cos(theta) ]   [ -T/m  ]   
            %
            %   omega_b x v_b = [q*zb_d - r*yb_d]
            %                   [r*xb_d - p*zb_d]
            %                   [p*yb_d - q*xb_d]
            %
            %   dv_b = F_b/m - (omega x v_b)
            %     xb_dd = r*yb_d - q*zb_d - g*sin(theta)          
            %     yb_dd = p*zb_d - r*xb_d + g*sin(phi)*cos(theta) 
            %     zb_dd = q*xb_d - p*yb_d + g*cos(phi)*cos(theta) - T/m 
            
            dxb_d = r*yb_d - q*zb_d - g*sthe;                  
            dyb_d = p*zb_d - r*xb_d + g*sphi*cthe;             
            dzb_d = q*xb_d - p*yb_d + g*cphi*cthe - (T/m);   


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
                    dxb_d; dyb_d; dzb_d;
                    dphi; dtheta; dpsi;
                    dp_; dq_; dr_];
        end


        function [A, B] = get_jacobians(obj, x_lin, u_lin)
            % GET_JACOBIANS  Analytical linearization at (x_lin, u_lin).
            %
            %   Returns A = df/dx (12x12) and B = df/du (12x4).
            %
            %
            %   At hover (phi=theta=psi=0, p=q=r=0, T=mg, v_b=0):
            %     A(1,4)  =  1         (Ẋ_I / xb_dot, R(1,1)=1)
            %     A(2,5)  =  1         (Ẏ_I / yb_dot)
            %     A(3,6)  =  1         (Ż_I / zb_dot)
            %     A(4,8)  = -g         (xb_dd / theta: gravity coupling)
            %     A(5,7)  = +g         (yb_dd / phi:   gravity coupling)
            %     B(6,1)  = -1/m
            %     B(10,2) = 1/Jxx   B(11,3)=1/Jyy   B(12,4)=1/Jzz

            % unpack
            xb_d  = x_lin(4);
            yb_d  = x_lin(5);
            zb_d  = x_lin(6);
            phi   = x_lin(7);
            theta = x_lin(8);
            psi   = x_lin(9);
            p     = x_lin(10);
            q     = x_lin(11);
            r     = x_lin(12);
            T     = u_lin(1);

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
            %% BLOCK 1 rows (1-3): P_I_dot = R(phi,theta,psi) * v_b
            %
            % Jacobian has TWO contributions:
            %   (a) d(R*v_b)/d(v_b) = R  [coupling body vel to inertial vel dot]
            %   (b) d(R*v_b)/d(angles) = dR/d(angle) * v_b  [coupling angles]
            %
            % Part (a): columns 4,5,6 get the rows of R
            %
            %   R row 1: [cthe*cpsi,  sphi*sthe*cpsi-cphi*spsi,  cphi*sthe*cpsi+sphi*spsi]
            %   R row 2: [cthe*spsi,  sphi*sthe*spsi+cphi*cpsi,  cphi*sthe*spsi-sphi*cpsi]
            %   R row 3: [-sthe,      sphi*cthe,                  cphi*cthe              ]

            % row 1 (dX_I), cols 4,5,6
            A(1,4) = cthe*cpsi;
            A(1,5) = sphi*sthe*cpsi - cphi*spsi;
            A(1,6) = cphi*sthe*cpsi + sphi*spsi;

            % row 2 (dY_I), cols 4,5,6
            A(2,4) = cthe*spsi;
            A(2,5) = sphi*sthe*spsi + cphi*cpsi;
            A(2,6) = cphi*sthe*spsi - sphi*cpsi;

            % row 3 (dZ_I), cols 4,5,6
            A(3,4) = -sthe;
            A(3,5) =  sphi*cthe;
            A(3,6) =  cphi*cthe;

            % Part (b): d(R*v_b)/d(phi), d(R*v_b)/d(theta), d(R*v_b)/d(psi)
            % [nonzero only off-hover when v_b != 0]
            %
            % d(dX_I)/dphi  = (cphi*sthe*cpsi+sphi*spsi)*yb_d + (-sphi*sthe*cpsi+cphi*spsi)*zb_d
            % d(dX_I)/dtheta= (-sthe*cpsi)*xb_d + (sphi*cthe*cpsi)*yb_d + (cphi*cthe*cpsi)*zb_d
            % d(dX_I)/dpsi  = (-cthe*spsi)*xb_d + (-sphi*sthe*spsi-cphi*cpsi)*yb_d + (-cphi*sthe*spsi+sphi*cpsi)*zb_d

            A(1,7) = ( cphi*sthe*cpsi + sphi*spsi)*yb_d + (-sphi*sthe*cpsi + cphi*spsi)*zb_d;
            A(1,8) = -sthe*cpsi*xb_d + sphi*cthe*cpsi*yb_d + cphi*cthe*cpsi*zb_d;
            A(1,9) = -cthe*spsi*xb_d + (-sphi*sthe*spsi - cphi*cpsi)*yb_d + (-cphi*sthe*spsi + sphi*cpsi)*zb_d;

            % d(dY_I)/dphi  = (cphi*sthe*spsi-sphi*cpsi)*yb_d + (-sphi*sthe*spsi-cphi*cpsi)*zb_d
            % d(dY_I)/dtheta= (-sthe*spsi)*xb_d + sphi*cthe*spsi*yb_d + cphi*cthe*spsi*zb_d
            % d(dY_I)/dpsi  = (cthe*cpsi)*xb_d + (sphi*sthe*cpsi-cphi*spsi)*yb_d + (cphi*sthe*cpsi+sphi*spsi)*zb_d

            A(2,7) = ( cphi*sthe*spsi - sphi*cpsi)*yb_d + (-sphi*sthe*spsi - cphi*cpsi)*zb_d;
            A(2,8) = -sthe*spsi*xb_d + sphi*cthe*spsi*yb_d + cphi*cthe*spsi*zb_d;
            A(2,9) =  cthe*cpsi*xb_d + (sphi*sthe*cpsi - cphi*spsi)*yb_d + (cphi*sthe*cpsi + sphi*spsi)*zb_d;

            % d(dZ_I)/dphi  = cphi*cthe*yb_d + (-sphi*cthe)*zb_d
            % d(dZ_I)/dtheta= -cthe*xb_d + (-sphi*sthe)*yb_d + (-cphi*sthe)*zb_d
            % d(dZ_I)/dpsi  = 0  (no psi dependence in Z_I_dot)

            A(3,7) =  cphi*cthe*yb_d - sphi*cthe*zb_d;
            A(3,8) = -cthe*xb_d - sphi*sthe*yb_d - cphi*sthe*zb_d;
            % A(3,9) = 0 by inspection


            %% BLOCK 2 rows (4-6): translational dynamics in body
            %
            %   dxb_d = r*yb_d - q*zb_d - g*sthe          
            %   dyb_d = p*zb_d - r*xb_d + g*sphi*cthe     
            %   dzb_d = q*xb_d - p*yb_d + g*cphi*cthe - T/m 
            %
            % Derivatives w.r.t. body velocities:
            %   d(dxb_d)/d(xb_d) = 0          d(dxb_d)/d(yb_d) = r   d(dxb_d)/d(zb_d) = -q
            %   d(dyb_d)/d(xb_d) = -r            d(dyb_d)/d(yb_d) = 0  d(dyb_d)/d(zb_d) = p
            %   d(dzb_d)/d(xb_d) =  q            d(dzb_d)/d(yb_d) = -p d(dzb_d)/d(zb_d) = 0
            %
            % Derivatives w.r.t. angles (gravity projection):
            %   d(dxb_d)/d(theta) = -g*cthe
            %   d(dyb_d)/d(phi)   =  g*cphi*cthe
            %   d(dyb_d)/d(theta) = -g*sphi*sthe
            %   d(dzb_d)/d(phi)   = -g*sphi*cthe
            %   d(dzb_d)/d(theta) = -g*cphi*sthe
            %
            % Derivatives w.r.t. rates (Coriolis):
            %   d(dxb_d)/d(q) = -zb_d   d(dxb_d)/d(r) =  yb_d
            %   d(dyb_d)/d(p) =  zb_d   d(dyb_d)/d(r) = -xb_d
            %   d(dzb_d)/d(p) = -yb_d   d(dzb_d)/d(q) =  xb_d

            % A(4,4)  = 0;
            A(4,5)  =  r;
            A(4,6)  = -q;
            A(4,8)  = -g*cthe;      % d(-g*sthe)/dtheta = -g*cthe
            A(4,11) = -zb_d;        % Coriolis: d(r*yb_d-q*zb_d)/dq = -zb_d
            A(4,12) =  yb_d;        % Coriolis: d/dr =  yb_d

            % row 5: dyb_d
            A(5,4)  = -r;
            % A(5,5)  = 0;
            A(5,6)  =  p;
            A(5,7)  =  g*cphi*cthe; % d(g*sphi*cthe)/dphi
            A(5,8)  = -g*sphi*sthe; % d(g*sphi*cthe)/dtheta
            A(5,10) =  zb_d;        % Coriolis: d(p*zb_d)/dp
            A(5,12) = -xb_d;        % Coriolis: d(-r*xb_d)/dr

            % row 6: dzb_d
            A(6,4)  =  q;
            A(6,5)  = -p;
            % A(6,6)  = 0;
            A(6,7)  = -g*sphi*cthe; % d(g*cphi*cthe)/dphi
            A(6,8)  = -g*cphi*sthe; % d(g*cphi*cthe)/dtheta
            A(6,10) = -yb_d;        % Coriolis
            A(6,11) =  xb_d;        % Coriolis


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
            %  Only dzb_d depends on T: d(dzb_d)/dT = -1/m

            B(6,1)  = -1/m;
            B(10,2) = 1/Jxx;
            B(11,3) = 1/Jyy;
            B(12,4) = 1/Jzz;

        end
    end
end