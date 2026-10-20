% SHAPECIRCLEDEPLOYMENT  Circle that starts from a hover setpoint, with tangent yaw.
%
% The drone is parked at (x0, y0, z0) with heading psi0. It
% accelerates forward along its nose, turns right (psi increasing) on a circle of
% radius R, completes n_laps and stops again at the start point.
%
% ONE scalar drives everything: phi(t) = angle travelled on the circle.
%   psi(t) = psi0 + phi(t)                                (yaw, tangent to the path)
%   x(t)   = x0 + R*( sin(psi(t)) - sin(psi0) )           (North)
%   y(t)   = y0 - R*( cos(psi(t)) - cos(psi0) )           (East)
%   z(t)   = z0                                           (Down, constant)
% Velocity = R*phi_dot*[cos(psi); sin(psi)]  ->  always along the nose.
%
% phi(t): angular speed goes 0 -> omega with an S-curve (T_ramp), stays at omega,
% goes omega -> 0 (T_ramp). S-curve s(u) = 35u^4 - 84u^5 + 70u^6 - 20u^7
% (s and its first 3 derivatives are 0 at both ends) -> at t = 0 and t = T all
% derivatives of x, y, psi vanish: phi = theta = 0, p = q = r = 0, T_i = m*g/4.
%
% Parameters (struct, e.g. TrajectoryParameters.CircleDeployment):
%   p0     [m]     [x0; y0; z0] start setpoint (NED)
%   psi0   [rad]   initial heading
%   R      [m]     radius
%   omega  [rad/s] cruise angular speed (cruise speed = R*omega)
%   n_laps [-]     integer -> start point = end point
%   T_ramp [s]     duration of speed-up and of slow-down
%
% Duration: duration() = n_laps*2*pi/omega + T_ramp.

classdef ShapeCircleDeployment < TrajectoryBase

    methods

        function obj = ShapeCircleDeployment(params)
            obj@TrajectoryBase(params);
            assert(params.T_ramp > 0 && params.T_ramp <= params.n_laps*2*pi/params.omega, ...
                'Need 0 < T_ramp <= n_laps * lap time.');
        end

        function T = duration(obj)
            T = obj.params.n_laps * 2*pi / obj.params.omega + obj.params.T_ramp;
        end

        function [sigma, dsigma, ddsigma, dddsigma, ddddsigma] = get_flat_outputs(obj, t)

            P  = obj.params;
            R  = P.R;
            x0 = P.p0(1);  y0 = P.p0(2);  z0 = P.p0(3);

            %% Yaw = psi0 + angle travelled (tangent to the circle)
            [phi, w1, w2, w3, w4] = obj.angleLaw(t);   % w_k = k-th derivative of phi (= of psi)
            psi = P.psi0 + phi;
            c = cos(psi);  s = sin(psi);

            %% Position and its derivatives (closed form)
            x     = x0 + R*(s - sin(P.psi0));
            y     = y0 - R*(c - cos(P.psi0));

            xd    = R*c*w1;
            yd    = R*s*w1;

            xdd   = R*(c*w2 - s*w1^2);
            ydd   = R*(s*w2 + c*w1^2);

            xddd  = R*(c*(w3 - w1^3) - 3*s*w1*w2);
            yddd  = R*(s*(w3 - w1^3) + 3*c*w1*w2);

            A4    = w4 - 6*w1^2*w2;
            B4    = 3*w2^2 + 4*w1*w3 - w1^4;
            xdddd = R*(c*A4 - s*B4);
            ydddd = R*(s*A4 + c*B4);

            %% ASSEMBLY  [x; y; z; psi]
            sigma     = [x;     y;     z0; psi];
            dsigma    = [xd;    yd;    0;  w1];
            ddsigma   = [xdd;   ydd;   0;  w2];
            dddsigma  = [xddd;  yddd;  0;  w3];
            ddddsigma = [xdddd; ydddd; 0;  w4];
        end

        function [phi, w1, w2, w3, w4] = angleLaw(obj, t)
            % phi(t) (angle travelled, phi(0) = 0) and its first four time derivatives
            P   = obj.params;
            w   = P.omega;
            Tr  = P.T_ramp;
            T   = obj.duration();
            ang = P.n_laps * 2*pi;
            t   = min(max(t, 0), T);

            if t <= Tr                                   % speed up
                [I, s0, s1, s2, s3] = ShapeCircleDeployment.scurve(t / Tr);
                phi = w*Tr*I;
                w1 = w*s0;  w2 = w*s1/Tr;  w3 = w*s2/Tr^2;  w4 = w*s3/Tr^3;
            elseif t >= T - Tr                           % slow down (mirror)
                [I, s0, s1, s2, s3] = ShapeCircleDeployment.scurve((T - t) / Tr);
                phi = ang - w*Tr*I;
                w1 = w*s0;  w2 = -w*s1/Tr;  w3 = w*s2/Tr^2;  w4 = -w*s3/Tr^3;
            else                                         % cruise
                phi = w*Tr/2 + w*(t - Tr);
                w1 = w;  w2 = 0;  w3 = 0;  w4 = 0;
            end
        end
    end

    methods (Static, Access = private)
        function [I, s0, s1, s2, s3] = scurve(u)
            % s0: S-curve 0 -> 1, s1..s3: derivatives w.r.t. u, I: integral (I(1) = 1/2)
            I  = 7*u^5 - 14*u^6 + 10*u^7 - 2.5*u^8;
            s0 = 35*u^4 - 84*u^5 + 70*u^6 - 20*u^7;
            s1 = 140*u^3 - 420*u^4 + 420*u^5 - 140*u^6;
            s2 = 420*u^2 - 1680*u^3 + 2100*u^4 - 840*u^5;
            s3 = 840*u - 5040*u^2 + 8400*u^3 - 4200*u^4;
        end
    end
end