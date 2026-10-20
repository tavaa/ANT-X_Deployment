% SHAPECIRCLE  Circular trajectory.
%
% Flat outputs:
%   sigma = [x_I; y_I; z_I; psi]
%
% Trajectory properties:
%   - Circular motion in the XY plane
%   - Constant altitude
%   - Constant angular speed
%   - Heading aligned with velocity direction

classdef ShapeCircle < TrajectoryBase

    methods

        function obj = ShapeCircle(circle_params)
            obj@TrajectoryBase(circle_params);
        end

        function [sigma, dsigma, ddsigma, dddsigma, ddddsigma] = ...
            get_flat_outputs(obj, t)

            %% Parameters
            R = obj.params.R;
            w = obj.params.omega;
            h = obj.params.z_height;

            wt = w * t;

            %% Position
            x = R * cos(wt);
            y = R * sin(wt);
            z = h;

            %% Velocity
            xd = -R * w * sin(wt);
            yd =  R * w * cos(wt);
            zd = 0;

            %% Accelleration
            xdd = -R * w^2 * cos(wt);
            ydd = -R * w^2 * sin(wt);
            zdd = 0;

            %% Jerk
            xddd =  R * w^3 * sin(wt);
            yddd = -R * w^3 * cos(wt);
            zddd = 0;

            %% Snap
            xdddd =  R * w^4 * cos(wt);
            ydddd =  R * w^4 * sin(wt);
            zdddd = 0;

            %% YAW
            % Heading aligned with velocity direction
            psi = wt + pi/2;
            psi_d   = w;
            psi_dd  = 0;
            psi_ddd = 0;
            psi_dddd = 0;

        
            %% ASSEMBLY
            sigma = [x; y; z; psi];
            dsigma = [xd; yd; zd; psi_d];
            ddsigma = [xdd; ydd; zdd; psi_dd];
            dddsigma = [xddd; yddd; zddd; psi_ddd];
            ddddsigma = [xdddd; ydddd; zdddd; psi_dddd];

        end
    end
end