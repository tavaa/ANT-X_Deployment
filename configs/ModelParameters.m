classdef ModelParameters
    % MODELPARAMETERS  Physical parameters of the ANT-X quadcopter.
    %
    % Sources:
    %   - v3: OFFICIAL lab parameters (parameters_drone.txt, tutor code).
    %   - v2 (commented): mass/arm from ANT-X spec sheet, inertia from
    %     Cavagnini (2021) thesis, C_T/C_Q by analogy, damping from
    %     Ghignoni closed-form estimate.
    %
    % v3 NOTE: aerodynamic damping/drag deliberately excluded (Lp=Mq=Nr=0).
    % Simple model sufficient for lab deployment.

properties (Constant)

        %% PHYSICAL PARAMETERS
        % m   = 0.270;   % [kg]    Total mass v1
        % m   = 0.363;   % [kg]    Total mass ANT-X v2 spec sheet
        m   = 0.291;     % [kg]    Total mass, OFFICIAL (lab) - v3
        b   = 0.08;      % [m]     Arm length, CoG to motor 
        g   = 9.81;      % [m/s^2] Gravitational acceleration

        rho = 1.225;     % [kg/m^3] Air density, sea level standard 
        R   = 0.0381;    % [m]     Propeller radius 
        D   = 0.0762;    % [m]     Propeller diameter 

        %% INERTIA TENSOR (Diagonal, X-config symmetry: Jxx = Jyy > Jzz)
        % Jxx = 0.00307;   % [kg*m^2] Roll inertia, Cavagnini (2021) v2
        % Jyy = 0.00307;   % [kg*m^2] Pitch inertia, Cavagnini (2021) v2
        % Jzz = 0.00239;   % [kg*m^2] Yaw inertia, Cavagnini (2021) v2
        Jxx = 0.00067;     % [kg*m^2] Roll inertia, OFFICIAL (lab) v3
        Jyy = 0.00067;     % [kg*m^2] Pitch inertia, OFFICIAL (lab) v3
        Jzz = 0.00116;     % [kg*m^2] Yaw inertia, OFFICIAL (lab) v3

        %% PROPULSION CONSTANTS
        % C_T = 0.0343;    % [-] Thrust coefficient, preliminary v2
        % C_Q = 0.00841;   % [-] Torque coefficient, preliminary v2
        % K_T  = 2.7823e-07;   % [N*s^2/rad^2]   Thrust constant, v2
        % K_Q  = 2.5994e-09;   % [N*m*s^2/rad^2] Torque constant, v2
        % Beta = 0.009340;     % [m] Torque-to-thrust ratio, v2

        K_T  = 3.168e-07;      % [N*s^2/rad^2]   Thrust constant v3
        Beta = 0.01067;        % [m] Torque-to-thrust ratio v3
        K_Q  = 3.380e-09;   % [N*m*s^2/rad^2] Torque constant v3
        Omega_max = 2680;      % [rad/s] Max rotor speed v3

        %% ROTATIONAL DAMPING
        % Lp = -7.31e-04;   % [Nm*s] Ghignoni closed-form estimate, preliminary
        % Mq = -7.31e-04;   % [Nm*s] Ghignoni closed-form estimate, preliminary
        Lp = 0;            % [Nm*s] excluded in v3
        Mq = 0;            % [Nm*s] excluded in v3 
        Nr = 0;            % [Nm*s] excluded in v3

        %% TRANSLATIONAL DRAG - ignored
        % X_u = ;
        % Y_v = ;
        % Z_w = ;

        %% MOTOR CONSTANTS
        % tau_accel = 0.035;   % [s] Accel-phase time constant, v2
        % tau_decel = 0.035;   % [s] Decel-phase time constant, v2
        tau_accel = 0.0125;    % [s] Accel-phase time constant v3
        tau_decel = 0.025;     % [s] Decel-phase time constant v3
        tau_p = 0.01875;       % [s] Motor-lag time constant (symmetric) tau_p = ((tau_acel + tau_decel)/2)
end
end