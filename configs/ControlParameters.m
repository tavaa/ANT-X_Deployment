% CONTROLPARAMS Parameters defines the Controllers and system's sampling
% time properties
% Defines the State and Inputs constraints, optimization horizon, sampling rate, and weighting matrices.
%
% Q_diag >= 0 (Positive semi-definite state cost).
% R_diag > 0 (Positive definite input cost for solver stability).

classdef ControlParameters
    
    properties (Constant)

        % Sampling time [s] Configuration
        Ts = 0.10;    % - 10Hz

        %% State and Inputs Constraints
        % State Constraints: Bounds for both models
        % A00_SimplifiedModel - A02_SimplifiedModel_Deployment
        % TODO modify state constraints here (inertal + check others)
        x_min_12 = [-inf; -inf; -inf; -3; -3; -3; -deg2rad(35); -deg2rad(35); -inf; -deg2rad(600); -deg2rad(600); -deg2rad(600)];
        x_max_12 = [inf; inf; inf; 3; 3; 1; +deg2rad(35); +deg2rad(35); inf; deg2rad(600); deg2rad(600); deg2rad(600)];
        
        % A01_SimplifiedModel_quat
        x_min_13 = [-inf; -inf; -inf; -3; -3; -3; -inf; -inf; -inf; -inf; -deg2rad(600); -deg2rad(600); -deg2rad(600)];
        x_max_13 = [inf; inf; inf; 3; 3; 1; +inf; +inf; +inf; +inf; deg2rad(600); deg2rad(600); deg2rad(600)];

        % B01_ANTX_quadcopter
        x_min_16 = [-inf; -inf; -inf; -3; -3; -3; -deg2rad(35); -deg2rad(35); -inf; -deg2rad(600); -deg2rad(600); -deg2rad(600); 0; 0; 0; 0];
        x_max_16 = [inf; inf; inf; 3; 3; 1; +deg2rad(35); +deg2rad(35); inf; deg2rad(600); deg2rad(600); deg2rad(600); 2.27538; 2.27538; 2.27538; 2.27538];
        
        %% V1 - Constraints on virtual inputs
        % Input Saturation Limits [Thrust, Roll, Pitch, Yaw]
        % u in [u_min, u_max].
        u_min = [0.0; -0.15; -0.15; -0.05]; 
        u_max = [ 2 * 0.270 * 9.81;  0.15;  0.15;  0.05]; 

        % deployment input constraints:
        % acts on [T, L, M, N]
        u_min_12 = [0.0; 0.0; 0.0; 0.0] 
        u_max_12 = [9.10; 0.257; 0.257; 0.0486];

        %% V3 - Constraints on Omega [rad/s] for each rotor
        % u in R^4 = desired motor thrusts T_i^d [N].
        %   Omega_hover = sqrt(m*g/(4*K_T))
        %   Omega_max   = 2680 rad/s
        %   Omega_min   = 0,    
        %   T_max = K_T*Omega_max^2,  T_min = K_T*Omega_min^2
        u_min_16 = [0; 0; 0; 0];
        u_max_16 = [2.27538; 2.27538; 2.27538; 2.27538];
        
        %% MPC Parameters
        % Prediction horizon
        p  = 18;  
        
        % State Weighting Matrix Q (Diagonal entries)
        % Q_diag_12 -> 12 states models - Euler's angles
        % Q_diag_13 -> 13 states models - Quaternions
        % Vector 12 states: [x_I, y_I, z_I, vx_B, vy_B, vz_B, phi, theta, psi, p, q, r] 
        % Vector 13 states: [x_I, y_I, z_I, vx_B, vy_B, vz_B, q0, q1, q2, q3, p, q, r] 
        % Vector 16 states: [r_n, r_e, r_d, v_d, v_e, v_d, phi, theta, psi, p, q, r, T_1, T_2, T_3, T_4] 
        Q_diag_12 = [50; 50; 50; 8; 8; 12; 8; 8; 10; 1; 1; 1];
        Qf_diag_12 = [50; 50; 50; 8; 8; 12; 8; 8; 10; 1; 1; 1]; % A02 -> Q_f = 10 * Q;

        Q_diag_13  = [50; 50; 50; 8; 8; 12; 0; 32; 32; 40; 1; 1; 1];
        Qf_diag_13 = [50; 50; 50; 8; 8; 12; 0; 32; 32; 40; 1; 1; 1];

        % V2 -> tuned using Bryson Rule on state/input constraints
        %Q_diag_16  = [50; 50; 50; 0.11; 0.11; 1; 2.67; 2.67; 10; 0.1; 0.1; 0.1; 2.0; 2.0; 2.0; 2.0];
        %Qf_diag_16  = 50 * [50; 50; 50; 0.11; 0.11; 1; 2.67; 2.67; 10; 0.1; 0.1; 0.1; 2.0; 2.0; 2.0; 2.0];

        % Input Weighting Matrix R (Diagonal entries)
        % Vector: [T, L, M, N]
        R_diag_12 = [2.0; 2.0; 2.0; 2.0];
        % Vector: [T_1^d, T_2^d, T_3^d, T_4^d]
        %R_diag_16 = [2.0; 2.0; 2.0; 2.0];
        
        % V3 -> tuned using Bryson Rule on desired state/input error (more aggressive, to test in lab.) 
        Q_diag_16  = [100;100;100; 25;25;25; 131.31;131.31;131.31; 4;4;4; 400;400;400;400];
        Qf_diag_16 = 10 * [100;100;100; 25;25;25; 131.31;131.31;131.31; 4;4;4; 400;400;400;400];
        R_diag_16  = [100; 100; 100; 100];
    end
end

