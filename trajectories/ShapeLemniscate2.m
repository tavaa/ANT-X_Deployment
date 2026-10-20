classdef ShapeLemniscate2 < TrajectoryBase
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

        function obj = ShapeLemniscate2(params)
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

            wt  = w * t;
            wt2 = 2 * w * t;

            %% TRANSLATION (Gerono Lemniscate)
            % Parameterization:
            %   x(t) = R * cos(w*t)
            %   y(t) = R * sin(w*t)*cos(w*t) = (R/2) * sin(2*w*t)
            
            % Position
            x = R * cos(wt) + cx;
            y = 0.5 * R * sin(wt2) + cy;
            z = h;

            % Velocity (1st Derivative)
            xd = -R * w * sin(wt);
            yd =  R * w * cos(wt2);
            zd = 0;

            % Acceleration (2nd Derivative)
            xdd = -R * w^2 * cos(wt);
            ydd = -2 * R * w^2 * sin(wt2);
            zdd = 0;

            % Jerk (3rd Derivative)
            xddd =  R * w^3 * sin(wt);
            yddd = -4 * R * w^3 * cos(wt2);
            zddd = 0;

            % Snap (4th Derivative)
            xdddd = R * w^4 * cos(wt);
            ydddd = 8 * R * w^4 * sin(wt2);
            zdddd = 0;


            %% YAW

            %Heading aligned with velocity vector
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