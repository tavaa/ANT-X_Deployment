classdef ShapeYawTest < TrajectoryBase
    % SHAPEYAWTEST  Pure yaw at fixed position.
    %
    % Flat outputs:
    %   sigma = [x_I; y_I; z_I; psi]
    %
    % Trajectory properties:
    %   - Hover at fixed inertial position
    %   - Smooth yaw-only excitation
    %   - Isolates yaw rotational dynamics
    %   - C^4 continuous yaw profile (snap continuous)
    %
    % Parameters
    %   A        [rad]    Yaw amplitude
    %   omega    [rad/s]  Trajectory frequency
    %   z_height [m NED]  Constant altitude

    methods

        function obj = ShapeYawTest(params)
            obj@TrajectoryBase(params);
        end

        function [sigma, dsigma, ddsigma, dddsigma, ddddsigma] = ...
                get_flat_outputs(obj, t)
            %GET_FLAT_OUTPUTS Generate smooth yaw-only flat outputs.

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

            %% Fixed position
            x = 0;
            y = 0;

            %% Zero translational derivatives
            dx = 0;
            dy = 0;
            dz = 0;

            ddx = 0;
            ddy = 0;
            ddz = 0;

            dddx = 0;
            dddy = 0;
            dddz = 0;

            ddddx = 0;
            ddddy = 0;
            ddddz = 0;

            %% Smooth yaw trajectory
            % psi(t) = K * (1 - cos(wt))^3
            psi = K * (1 - c)^3;

            %% Yaw rate
            dpsi = 3 * K * w * s * (1 - c)^2;

            %% Yaw acceleration
            ddpsi = 3 * K * om2 * (1 - c)^2 * (3 * c + 2);

            %% Yaw jerk
            dddpsi = 3 * K * om3 * s * (1 - c) * (9 * c + 1);

            %% Yaw snap
            ddddpsi = 3 * K * om4 * (1 - c) * (27 * c^2 + 11 * c - 8);

            %% Pack flat outputs
            sigma = [x; y; z; psi];
            dsigma = [dx; dy; dz; dpsi];
            ddsigma = [ddx; ddy; ddz; ddpsi];
            dddsigma = [dddx; dddy; dddz; dddpsi];
            ddddsigma = [ddddx; ddddy; ddddz; ddddpsi];

        end
    end
end