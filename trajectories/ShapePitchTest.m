classdef ShapePitchTest < TrajectoryBase

    methods

        function obj = ShapePitchTest(params)
            obj@TrajectoryBase(params);
        end

        function [sigma, dsigma, ddsigma, dddsigma, ddddsigma] = ...
                get_flat_outputs(obj, t)
            %GET_FLAT_OUTPUTS Generate smooth longitudinal flat-output trajectory.
            %
            % Flat outputs:
            %   sigma = [x_I; y_I; z_I; psi]
            %
            % Trajectory properties:
            %   - Motion only along inertial x-axis
            %   - y and yaw remain constant
            %   - z fixed at constant altitude
            %   - C^4 smooth profile (snap continuous)

            %% Parameters
            A = obj.params.A;
            w = obj.params.omega;
            z = obj.params.z_height;

            %% Common terms
            K = A / 8;

            wt = w * t;
            c = cos(wt);
            s = sin(wt);

            om2 = w^2;
            om3 = w^3;
            om4 = w^4;

            %% Position trajectory (x-axis)
            %
            % Smooth shaping polynomial:
            %
            %   x(t) = K * (1 - cos(wt))^3

            x = K * (1 - c)^3;

            %% Velocity
            dx = 3 * K * w * s * (1 - c)^2;

            %% Acceleration
            ddx = 3 * K * om2 * (1 - c)^2 * (3 * c + 2);

            %% Jerk
            dddx = 3 * K * om3 * s * (1 - c) * (9 * c + 1);

            %% Snap
            ddddx = 3 * K * om4 * (1 - c) * (27 * c^2 + 11 * c - 8);

            %% Constant axes
            y     = 0;
            dy    = 0;
            ddy   = 0;
            dddy  = 0;
            ddddy = 0;

            %% Constant heading
            psi     = 0;
            dpsi    = 0;
            ddpsi   = 0;
            dddpsi  = 0;
            ddddpsi = 0;

            %% Pack flat outputs
            sigma = [x; y; z; psi];
            dsigma = [dx; dy; 0; dpsi];
            ddsigma = [ddx; ddy; 0; ddpsi];
            dddsigma = [dddx; dddy; 0; dddpsi];
            ddddsigma = [ddddx; ddddy; 0; ddddpsi];

        end
    end
end