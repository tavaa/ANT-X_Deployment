classdef B01_Allocator
    % B01_ALLOCATOR  Virtual/physical allocation utilities for the ANT-X
    %                quadrotor, B01 motor convention.
    %
    % Virtual inputs:  [T; L; M; N]        (collective thrust, roll/pitch/yaw moments)
    % Physical inputs: [T1; T2; T3; T4]    (individual motor thrusts)
    %
    %   virtual  = Chi * physical
    %   physical = Chi \ virtual

    properties
        mp   % ModelParameters 
    end

    methods
        function obj = B01_Allocator(mp)
            arguments
                mp (1,1) ModelParameters
            end
            obj.mp = mp;
        end

        function C = Chi(obj)
            % Mixer matrix: [T;L;M;N] = C * [T1;...;T4]
            b    = obj.mp.b;
            beta = obj.mp.Beta;
            C = [ 1,          1,          1,          1;
                  b/sqrt(2),  b/sqrt(2), -b/sqrt(2), -b/sqrt(2);
                  b/sqrt(2), -b/sqrt(2), -b/sqrt(2),  b/sqrt(2);
                  beta,      -beta,       beta,      -beta ];
        end

        function T_motors = virtualToThrust(obj, virtual)
            % [T;L;M;N] -> [T1;...;T4]
            T_motors = obj.Chi() \ virtual;
        end

        function omega_sq = virtualToOmegaSq(obj, virtual)
            % [T;L;M;N] -> [Omega1^2;...;Omega4^2]
            omega_sq = obj.virtualToThrust(virtual) / obj.mp.K_T;
        end

        function omega = virtualToOmega(obj, virtual)
            % [T;L;M;N] -> [Omega1;...;Omega4], negative Omega^2 clipped to 0
            omega_sq = obj.virtualToOmegaSq(virtual);
            omega_sq(omega_sq < 0) = 0;
            omega = sqrt(omega_sq);
        end

        function virtual = omegaSqToVirtual(obj, omega_sq)
            % [Omega1^2;...;Omega4^2] -> [T;L;M;N]
            virtual = obj.Chi() * (obj.mp.K_T * omega_sq);
        end
    end
end