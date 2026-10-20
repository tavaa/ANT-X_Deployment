% TRAJECTORYBASE Abstract Interface for Differentially Flat Trajectory Generation
%
% Flat outputs:
%   sigma = [x_I; y_I; z_I; psi]
%
% Required derivatives:
%   dsigma      -> velocity / yaw rate
%   ddsigma     -> acceleration / yaw acceleration
%   dddsigma    -> jerk / yaw jerk
%   ddddsigma   -> snap / yaw snap

classdef (Abstract) TrajectoryBase < handle

    properties
        % params : trajectory parameter struct
        params
    end

    methods

        function obj = TrajectoryBase(params)
            obj.params = params;
        end

    end

    methods (Abstract)

        % GET_FLAT_OUTPUTS
        %
        % Returns flat outputs and derivatives up to 4th order.
        %
        % Inputs
        %   t : time [s]
        %
        % Outputs
        %   sigma       : [4x1]
        %   dsigma      : [4x1]
        %   ddsigma     : [4x1]
        %   dddsigma    : [4x1]
        %   ddddsigma   : [4x1]

        [sigma, dsigma, ddsigma, dddsigma, ddddsigma] = ...
            get_flat_outputs(obj, t);

    end

end