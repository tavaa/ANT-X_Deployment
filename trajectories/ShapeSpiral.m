% SHAPESPIRAL Spiral trajectory.
%
% Flat outputs:
%   sigma = [x_I; y_I; z_I; psi]
%
% Parameters 
%   R_x      [m]       Radius X
%   R_y      [m]       Radius Y
%   omega    [rad/s]   Angular Velocity
%   z_start  [m NED]   Initial Height
%   v_z      [m/s]     Climbing velocity (v_z > 0 = up, dz/dt = -v_z in NED)

classdef ShapeSpiral < TrajectoryBase

    methods

        function obj = ShapeSpiral(params)
            obj@TrajectoryBase(params);
        end

        function [sigma, dsigma, ddsigma, dddsigma, ddddsigma] = ...
                get_flat_outputs(obj, t)

            %% Parameters
            Rx = obj.params.Rx;
            Ry = obj.params.Ry;
            w  = obj.params.omega;
            z0 = obj.params.z_start;
            vz = obj.params.vz;

            wt = w * t;

            %% Position
            x = Rx * cos(wt);
            y = Ry * sin(wt);
            z = z0 - vz * t;

            %% Velocity
            xd = -Rx * w * sin(wt);
            yd =  Ry * w * cos(wt);
            zd = -vz;

            %% Accelleration
            xdd = -Rx * w^2 * cos(wt);
            ydd = -Ry * w^2 * sin(wt);
            zdd =  0;

            %% Jerk
            xddd =  Rx * w^3 * sin(wt);
            yddd = -Ry * w^3 * cos(wt);
            zddd =  0;

            %% Snap
            xdddd =  Rx * w^4 * cos(wt);
            ydddd =  Ry * w^4 * sin(wt);
            zdddd =  0;

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