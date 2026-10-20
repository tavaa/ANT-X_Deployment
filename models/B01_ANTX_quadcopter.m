classdef B01_ANTX_quadcopter
    % B01_ANTX_quadcopter  Nonlinear quadrotor model, X configuration.
    %
    % Newton-Euler physics, velocities expressed directly in the
    % INERTIAL NED frame. Includes rotational aerodynamic damping
    % (Lp, Mq, Nr). No translational drag.
    %
    % delay_motors FLAG (constructor argument, default true):
    %   true  -> 16-state model: adds 4 auxiliary states for first-order
    %            asymmetric motor delay. Input u = [T1_d;...;T4_d]
    %            (desired thrusts).
    %   false -> 12-state rigid-body-only model, no delay states.
    %            Input u = [T1;...;T4] directly (instantaneous thrust,
    %            T_i = T_i^d by construction).
    %
    % Reference frames
    %   Inertial : NED  (North-East-Down)
    %   Body     : FRD  (Forward-Right-Down)
    %
    % State vector, delay_motors = true, x in R^16 (column)
    %   x(1)  = r_n     [m]      North position     (NED)
    %   x(2)  = r_e     [m]      East  position     (NED)
    %   x(3)  = r_d     [m]      Down  position     (NED, positive downward)
    %   x(4)  = v_n     [m/s]    North velocity     (NED, inertial)
    %   x(5)  = v_e     [m/s]    East  velocity     (NED, inertial)
    %   x(6)  = v_d     [m/s]    Down  velocity     (NED, inertial)
    %   x(7)  = phi     [rad]    Roll  angle (Z-Y-X Euler)
    %   x(8)  = theta   [rad]    Pitch angle (positive nose-up)
    %   x(9)  = psi     [rad]    Yaw   angle
    %   x(10) = p       [rad/s]  Roll  rate (body)
    %   x(11) = q       [rad/s]  Pitch rate (body)
    %   x(12) = r       [rad/s]  Yaw   rate (body)
    %   x(13) = T1      [N]      Motor 1 actual thrust (Front-Left,  CCW)
    %   x(14) = T2      [N]      Motor 2 actual thrust (Rear-Left,   CW)
    %   x(15) = T3      [N]      Motor 3 actual thrust (Rear-Right,  CCW)
    %   x(16) = T4      [N]      Motor 4 actual thrust (Front-Right, CW)
    %   u     in R^4 = [T1_d; T2_d; T3_d; T4_d]   desired motor thrusts
    %
    % State vector, delay_motors = false, x in R^12 (column)
    %   x(1:12) same as above (rigid body only, no motor states)
    %   u       in R^4 = [T1; T2; T3; T4]   instantaneous motor thrusts
    %
    % Virtual inputs (collective thrust, roll/pitch/yaw moments) are NOT
    % states or inputs in either mode -- algebraic functions of T1..T4
    % (from x(13:16) if delay_motors, from u otherwise):
    %
    %   T_tot  = T1 + T2 + T3 + T4
    %   L_virt = (b/sqrt(2)) * ( T1 + T2 - T3 - T4)
    %   M_virt = (b/sqrt(2)) * ( T1 - T2 - T3 + T4)
    %   N_virt =  beta       * ( T1 - T2 + T3 - T4),   beta := K_Q/K_T
    %
    % Motor delay (delay_motors=true only): first-order ASYMMETRIC lag
    % on individual thrusts T1..T4 (not on Omega_i): tau_accel when
    % spinning up, tau_decel when spinning down, switched on
    % sign(Ti_d - Ti).
    %
    % Hover trim:
    %   delay_motors=true : x*=[zeros(12,1); (m*g/4)*ones(4,1)], u*=(m*g/4)*ones(4,1)
    %   delay_motors=false: x*=zeros(12,1), u*=(m*g/4)*ones(4,1)

    properties
        mp
        delay_motors (1,1) logical = true
    end

    methods

        function obj = B01_ANTX_quadcopter(params, delay_motors)
            % Constructor.
            %   params       : ModelParameters object (m, g,
            %                  Jxx, Jyy, Jzz, b, Beta, Lp, Mq, Nr, and,
            %                  if delay_motors, tau_accel, tau_decel)
            %   delay_motors : (optional, default true)
            arguments
                params
                delay_motors (1,1) logical = true
            end
            obj.mp = params;
            obj.delay_motors = delay_motors;
        end

        function dxdt = dynamics(obj, t, x, u)
            % DYNAMICS  Evaluate the nonlinear ODE dxdt = f(x,u).
            %
            %   Inputs
            %     t : scalar time [s]
            %     x : 16x1 (delay_motors=true) or 12x1 (false) state
            %     u : 4x1 input -- [T1_d;...;T4_d] if delay_motors,
            %                       [T1;...;T4]    otherwise
            %   Output
            %     dxdt : same size as x

            vn    = x(4);
            ve    = x(5);
            vd    = x(6);
            phi   = x(7);
            theta = x(8);
            psi   = x(9);
            p     = x(10);
            q     = x(11);
            r     = x(12);

            m   = obj.mp.m;
            g   = obj.mp.g;
            Jxx = obj.mp.Jxx;
            Jyy = obj.mp.Jyy;
            Jzz = obj.mp.Jzz;
            b     = obj.mp.b;
            beta  = obj.mp.Beta;
            Lp    = obj.mp.Lp;
            Mq    = obj.mp.Mq;
            Nr    = obj.mp.Nr;

            sphi = sin(phi);   cphi = cos(phi);
            sthe = sin(theta); cthe = cos(theta);
            spsi = sin(psi);   cpsi = cos(psi);

            assert(abs(cthe) > 1e-6, ...
                'QuadrotorModel:GimbalLock', ...
                'theta = %.2f rad: singular (gimbal lock).', theta);

            %% BLOCK 0 -- Virtual inputs from motor thrusts
            %
            % T1..T4 come from the STATE if delay_motors (actual, lagged
            % thrust), or directly from the INPUT if not (instantaneous).

            if obj.delay_motors
                T1 = x(13); T2 = x(14); T3 = x(15); T4 = x(16);
            else
                T1 = u(1);  T2 = u(2);  T3 = u(3);  T4 = u(4);
            end

            T_tot  = T1 + T2 + T3 + T4;
            L_virt = (b/sqrt(2)) * ( T1 + T2 - T3 - T4);
            M_virt = (b/sqrt(2)) * ( T1 - T2 - T3 + T4);
            N_virt =  beta       * ( T1 - T2 + T3 - T4);


            %% BLOCK 1 -- Position kinematics
            %
            % Velocities already expressed in NED, so trivial.

            dr_n = vn;
            dr_e = ve;
            dr_d = vd;


            %% BLOCK 2 -- Translational dynamics (Newton, inertial frame)
            %
            %   m*v_dot = gravity + R(phi,theta,psi) * [0;0;-T_tot]

            dvn = -(1/m) * (cpsi*sthe*cphi + spsi*sphi) * T_tot;
            dve = -(1/m) * (spsi*sthe*cphi - cpsi*sphi) * T_tot;
            dvd =  g - (1/m) * (cthe*cphi) * T_tot;


            %% BLOCK 3 -- Rotational kinematics (Euler Z-Y-X)

            tthe = sthe / cthe;
            iche = 1.0  / cthe;

            dphi   = p + sphi*tthe*q + cphi*tthe*r;
            dtheta =     cphi*q      - sphi*r;
            dpsi   =     sphi*iche*q + cphi*iche*r;


            %% BLOCK 4 -- Rotational dynamics with aerodynamic damping
            %
            %   Jxx*p_dot = (Jyy-Jzz)*q*r + L_virt + Lp*p
            %   Jyy*q_dot = (Jzz-Jxx)*r*p + M_virt + Mq*q
            %   Jzz*r_dot = (Jxx-Jyy)*p*q + N_virt + Nr*r

            dp_ = (1/Jxx) * ((Jyy - Jzz)*q*r + L_virt + Lp*p);
            dq_ = (1/Jyy) * ((Jzz - Jxx)*r*p + M_virt + Mq*q);
            dr_ = (1/Jzz) * ((Jxx - Jyy)*p*q + N_virt + Nr*r);

            dxdt_12 = [dr_n; dr_e; dr_d;
                       dvn; dve; dvd;
                       dphi; dtheta; dpsi;
                       dp_; dq_; dr_];

            if ~obj.delay_motors
                dxdt = dxdt_12;
                return;
            end


            %% BLOCK 5 -- Motor first-order delay (asymmetric, switched)
            %
            %   Ti_dot = -(1/tau_i)*Ti + (1/tau_i)*Ti_d
            %   tau_i = tau_accel  if Ti_d >= Ti  (spinning up)
            %         = tau_decel  if Ti_d <  Ti  (spinning down)

            tau_a = obj.mp.tau_accel;
            tau_d = obj.mp.tau_decel;

            T1d = u(1); T2d = u(2); T3d = u(3); T4d = u(4);

            if T1d >= T1, tau1 = tau_a; else, tau1 = tau_d; end
            if T2d >= T2, tau2 = tau_a; else, tau2 = tau_d; end
            if T3d >= T3, tau3 = tau_a; else, tau3 = tau_d; end
            if T4d >= T4, tau4 = tau_a; else, tau4 = tau_d; end

            dT1 = -(1/tau1)*T1 + (1/tau1)*T1d;
            dT2 = -(1/tau2)*T2 + (1/tau2)*T2d;
            dT3 = -(1/tau3)*T3 + (1/tau3)*T3d;
            dT4 = -(1/tau4)*T4 + (1/tau4)*T4d;

            dxdt = [dxdt_12; dT1; dT2; dT3; dT4];
        end


        function [A, B] = get_jacobians(obj, x_lin, u_lin)
            % GET_JACOBIANS  Analytical linearization at (x_lin, u_lin).
            %
            %   Returns A = df/dx, B = df/du.
            %   delay_motors=true  : A is 16x16, B is 16x4.
            %   delay_motors=false : A is 12x12, B is 12x4.

            phi   = x_lin(7);
            theta = x_lin(8);
            psi   = x_lin(9);
            p     = x_lin(10);
            q     = x_lin(11);
            r     = x_lin(12);

            m   = obj.mp.m;
            Jxx = obj.mp.Jxx;
            Jyy = obj.mp.Jyy;
            Jzz = obj.mp.Jzz;
            b     = obj.mp.b;
            beta  = obj.mp.Beta;
            Lp    = obj.mp.Lp;
            Mq    = obj.mp.Mq;
            Nr    = obj.mp.Nr;

            sphi = sin(phi);   cphi = cos(phi);
            sthe = sin(theta); cthe = cos(theta);
            spsi = sin(psi);   cpsi = cos(psi);
            tthe = sthe / cthe;
            iche = 1.0  / cthe;

            assert(abs(cthe) > 1e-6, ...
                'QuadrotorModel:GimbalLock', ...
                'theta = %.2f rad: singular (gimbal lock).', theta);

            % T1..T4 at the linearization point -- from state if
            % delay_motors, from input otherwise.
            if obj.delay_motors
                T1 = x_lin(13); T2 = x_lin(14); T3 = x_lin(15); T4 = x_lin(16);
            else
                T1 = u_lin(1);  T2 = u_lin(2);  T3 = u_lin(3);  T4 = u_lin(4);
            end
            T_tot = T1 + T2 + T3 + T4;

            n = 12 + 4*obj.delay_motors;
            A = zeros(n, n);
            B = zeros(n, 4);

            %% BLOCK 1 rows (1-3): position kinematics -- trivial

            A(1,4) = 1;
            A(2,5) = 1;
            A(3,6) = 1;


            %% BLOCK 2 rows (4-6): translational dynamics
            %
            % T1..T4 derivatives go to A columns 13-16 if delay_motors
            % (state), or to B columns 1-4 if not (input) -- same
            % closed-form value either way, just relocated.

            A(4,7) =  (T_tot/m) * (cpsi*sthe*sphi - spsi*cphi);
            A(4,8) = -(T_tot/m) * cpsi*cthe*cphi;
            A(4,9) =  (T_tot/m) * (spsi*sthe*cphi - cpsi*sphi);

            A(5,7) =  (T_tot/m) * (spsi*sthe*sphi + cpsi*cphi);
            A(5,8) = -(T_tot/m) * spsi*cthe*cphi;
            A(5,9) = -(T_tot/m) * (cpsi*sthe*cphi + spsi*sphi);

            A(6,7) = (T_tot/m) * cthe*sphi;
            A(6,8) = (T_tot/m) * sthe*cphi;
            % A(6,9) = 0  (no psi dependence)

            dv_dT = [ -(1/m) * (cpsi*sthe*cphi + spsi*sphi);
                      -(1/m) * (spsi*sthe*cphi - cpsi*sphi);
                      -(1/m) * cthe*cphi ];   % [d(dvn);d(dve);d(dvd)] / dTi, same for i=1..4

            if obj.delay_motors
                A(4,13:16) = dv_dT(1);
                A(5,13:16) = dv_dT(2);
                A(6,13:16) = dv_dT(3);
            else
                B(4,1:4) = dv_dT(1);
                B(5,1:4) = dv_dT(2);
                B(6,1:4) = dv_dT(3);
            end


            %% BLOCK 3 rows (7-9): Euler kinematics

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


            %% BLOCK 4 rows (10-12): rotational dynamics with damping
            %
            % L_virt/M_virt/N_virt derivatives w.r.t. T1..T4 go to A
            % columns 13-16 (state) or B columns 1-4 (input), same
            % closed-form value, relocated as above.

            A(10,11) = (1/Jxx) * (Jyy - Jzz) * r;
            A(10,12) = (1/Jxx) * (Jyy - Jzz) * q;
            A(10,10) = (1/Jxx) * Lp;

            A(11,10) = (1/Jyy) * (Jzz - Jxx) * r;
            A(11,12) = (1/Jyy) * (Jzz - Jxx) * p;
            A(11,11) = (1/Jyy) * Mq;

            A(12,10) = (1/Jzz) * (Jxx - Jyy) * q;
            A(12,11) = (1/Jzz) * (Jxx - Jyy) * p;
            A(12,12) = (1/Jzz) * Nr;

            dL_dT = (b/sqrt(2)) * [ 1,  1, -1, -1] / Jxx;
            dM_dT = (b/sqrt(2)) * [ 1, -1, -1,  1] / Jyy;
            dN_dT =  beta       * [ 1, -1,  1, -1] / Jzz;

            if obj.delay_motors
                A(10,13:16) = dL_dT;
                A(11,13:16) = dM_dT;
                A(12,13:16) = dN_dT;
            else
                B(10,1:4) = dL_dT;
                B(11,1:4) = dM_dT;
                B(12,1:4) = dN_dT;
            end

            if ~obj.delay_motors
                return;   % 12x12 A, 12x4 B -- done
            end


            %% BLOCK 5 rows (13-16): motor delay -- diagonal only

            tau_a = obj.mp.tau_accel;
            tau_d = obj.mp.tau_decel;
            T1d = u_lin(1); T2d = u_lin(2); T3d = u_lin(3); T4d = u_lin(4);

            if T1d >= T1, tau1 = tau_a; else, tau1 = tau_d; end
            if T2d >= T2, tau2 = tau_a; else, tau2 = tau_d; end
            if T3d >= T3, tau3 = tau_a; else, tau3 = tau_d; end
            if T4d >= T4, tau4 = tau_a; else, tau4 = tau_d; end

            A(13,13) = -1/tau1;
            A(14,14) = -1/tau2;
            A(15,15) = -1/tau3;
            A(16,16) = -1/tau4;

            B(13,1) = 1/tau1;
            B(14,2) = 1/tau2;
            B(15,3) = 1/tau3;
            B(16,4) = 1/tau4;

        end
    end
end