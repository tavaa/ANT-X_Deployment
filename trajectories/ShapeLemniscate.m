classdef ShapeLemniscate < TrajectoryBase
% SHAPELEMNISCATE trajectory.
%
% Flat outputs:
%   sigma = [x_I; y_I; z_I; psi]
%
% Parameters 
%   R        [m]       Radius X
%   omega    [rad/s]   Angular Velocity
%   z_height [m NED]   Initial Height
%   center   [m XY]    Position of the Lemniscate (intersection)

    methods

        function obj = ShapeLemniscate(params)
            if ~isfield(params, 'center')
                params.center = [0.0, 0.0];
            end
            obj@TrajectoryBase(params);
        end

        function [sigma, dsigma, ddsigma, dddsigma, ddddsigma] = ...
                get_flat_outputs(obj, t)

            %% Parameters
            R  = obj.params.R;
            w  = obj.params.omega;
            h  = obj.params.z_height;
            cx = obj.params.center(1);
            cy = obj.params.center(2);

            wt = w * t;
            A  = R * sqrt(2);  

            % trigonometric functions
            s  = sin(wt);
            c  = cos(wt);
            s2 = s^2;
            c2 = c^2;
            D  = 1 + s2;       

            %% Position
            x = A * c / D  + cx;
            y = A * s * c / D + cy;
            z = h;

            %% Velocity 
            xd = A*w * (s^3 - 3 * s) / D^2;
            yd =  A*w * (1 - 3*s2) / D^2;
            zd = 0;

            %% Accelleration
            xdd = A*w^2 * c * (-3 + 12 * s2 - s^4) / D^3;
            ydd = 2*A*w^2 * s * c * (3*s2 - 5) / D^3;
            zdd = 0;

            %% Jerk
            xddd = A * w^3 * s * (45 - 103 * s^2 + 43 * s^4 - s^6) / D^4;
            yddd = A * w^3 * (-10 + 88 * s^2 - 82 * s^4 + 12 * s^6) / D^4;
            zddd = 0;

            %% Snap
            xdddd = A * w^4 * c * (45 - 624 * s^2 + 730 * s^4 - 136 * s^6 + s^8) / D^5;
            ydddd = A * w^4 * 8 * s * c * (32 - 107 * s^2 + 50 * s^4 - 3 * s^6) / D^5;
            zdddd = 0;

            %% YAW
            % Heading aligned with velocity vector
            v2    = xd^2 + yd^2;
            eps_v = 1e-8;
            psi   = atan2(yd, xd);

             if v2 < eps_v
                psi_d    = 0;
                psi_dd   = 0;
                psi_ddd  = 0;
                psi_dddd = 0;

             else
                N1    = xd*ydd  - yd*xdd;
                psi_d = N1 / v2;
                v2_d   = 2*(xd*xdd + yd*ydd);
                N1_d   = xd*yddd - yd*xddd;
                
                psi_dd = (N1_d*v2 - N1*v2_d) / v2^2;
                psi_ddd = 0;
                psi_dddd = 0;

            end

            %% ASSEMBLY
            sigma = [x; y; z; psi];
            dsigma = [xd; yd; zd; psi_d];
            ddsigma = [xdd; ydd; zdd; psi_dd];
            dddsigma = [xddd; yddd; zddd; psi_ddd];
            ddddsigma = [xdddd; ydddd; zdddd; psi_dddd];

        end
    end
end