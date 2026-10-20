% TRAJECTORYPARAMETERS Config class for reference trajectory generation.
% All spatial coordinates use the NED frame (z < 0 is above ground).
classdef TrajectoryParameters
    
    properties (Constant)

        %% LAP DURATION [s]
        T_lap = 12.6;  

        %% HOVER
        Hover = struct( ...
            'pos', [0.0; 0.0; -1.0], ... % [m] NED pos (1m up)
            'psi', 0.0 ...               % [rad] Yaw angle
        );

        %% SINGLE DOF TESTS (Oscillatory)
        % Pitch test
        Pitch = struct( ...
              'A', 1.0, ...              % [m] Amplitude
              'omega', 1.0, ...          % [rad/s] Freq
              'z_height', -1.0 ...       % [m] NED alt
        );
        
        % Roll test
        Roll = struct( ...
              'A', 1.0, ...              % [m] Amplitude
              'omega', 1.0, ...          % [rad/s] Freq
              'z_height', -1.0 ...       % [m] NED alt
        );
        
        % Yaw test
        Yaw = struct( ...
              'A', +pi/2, ...             % [rad] Amplitude
              'omega', 1.0, ...          % [rad/s] Freq
              'z_height', -1.0 ...       % [m] NED alt
        );
        
        % Thrust test
        Thrust = struct( ...
              'H', 1.0, ...              % [m] Climb height
              'omega', 1.0, ...          % [rad/s] Freq
              'z_ground', -1.0 ...      % [m] NED start alt
        );

        %% 2D Shapes - XY 
        % Circular trajectory
        Circle = struct( ...
            'R', 1.0, ...                % [m] Radius
            'omega', 2*pi/TrajectoryParameters.T_lap, ...            % [rad/s] Angular speed
            'z_height', -2.0, ...        % [m] NED alt
            'center', [0.0; 0.0] ...     % [m] XY center
        );
        
        CircleDeployment = struct( ...
            'p0',     [0.0; 0.0; -1.0], ...   % [m] setpoint di partenza (NED)
            'psi0',   0.0, ...                % [rad] muso iniziale
            'R',      1.0, ...                % [m]
            'omega',  2*pi/TrajectoryParameters.T_lap, ...   % velocità di crociera R*omega = 0.5 m/s
            'n_laps', 1, ...
            'T_ramp', 3.0);                   % [s]

        % Lemniscate trajectory
        Lemniscate = struct( ...
            'R', 1.0, ...                % [m] Max amplitude scale
            'omega', 2*pi/TrajectoryParameters.T_lap, ...            % [rad/s] Angular speed
            'z_height', -2.0, ...        % [m] NED alt
            'center', [0.0; 0.0] ...     % [m] XY center
        );

        % Lemniscate2 trajectory
        Lemniscate2 = struct( ...
            'R', 1.0, ...                % [m] Max amplitude scale
            'omega', 2*pi/TrajectoryParameters.T_lap, ...            % [rad/s] Angular speed
            'z_height', -2.0, ...        % [m] NED alt
            'center', [0.0; 0.0] ...     % [m] XY center
        );

        %% 3D Shapes - XYZ
        % Spiral Trajectory
        Spiral = struct( ...
            'Rx', 1.5, ...               % [m] X-axis radius
            'Ry', 1.0, ...               % [m] Y-axis radius
            'omega', 2*pi/TrajectoryParameters.T_lap, ...            % [rad/s] Angular speed
            'vz', 0.05, ...             % [m/s] Climb rate
            'z_start', -1.5 ...          % [m] NED start alt
        );

        % 3D Lissajous trajectory 
        Lissajous3D = struct( ... 
            'Ax', 1.0, ... % [m] X amplitude 
            'Ay', 0.5, ... % [m] Y amplitude 
            'Az', 0.20, ... % [m] Z amplitude 
            'omega', 2*pi/TrajectoryParameters.T_lap, ... % [rad/s] Angular speed 
            'z_height', -2.0, ... % [m] Mean NED altitude
            'psi0', 0.0, ... % [rad] Initial heading 
            'center', [0.0; 0.0] ... % [m] XY center 
         );


        %% SIMULATION DURATIONS [s]
        % T = (n_laps * 2*pi) / omega
        T_duration_hover      = 5.0;            % Fixed duration
        T_duration_pitch      = (2*pi / 1.0);     % 1 period
        T_duration_roll       = (2*pi / 1.0);     % 1 period
        T_duration_yaw        = (2*pi / 1.0);     % 1 period
        T_duration_thrust     = (2*pi / 1.0);     % 1 period
        T_duration_circle      = 2 * TrajectoryParameters.T_lap;   % N_laps * circle lap
        T_duration_circle_deployment      = 2 * TrajectoryParameters.T_lap + 3.0;   % N_laps * circle lap
        T_duration_lemniscate  = 2 * TrajectoryParameters.T_lap;   % N_laps * figure-8 lap
        T_duration_lemniscate2 = 2 * TrajectoryParameters.T_lap;   % N_laps * figure-8 lap
        T_duration_spiral      = 2 * TrajectoryParameters.T_lap;   % N_laps * spiral lap
        T_duration_lissajous3D = 2 * TrajectoryParameters.T_lap;   % N_laps * lissajous lap
        
    end
end