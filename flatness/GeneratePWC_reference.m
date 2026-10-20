function [x_ref, u_ref, t_steps] = GeneratePWC_reference(traj_obj, T_sim, Ts, model, flatness, options, quad_options)
% GENERATEPWC_REFERENCE  Generates a Piece-Wise Constant (PWC) reference
%                        by integrating the nonlinear model forward with
%                        inputs held constant at each sampling instant.
%
% Supports both 12-state (Euler Z-Y-X) and 13-state (Quaternion) representations.
%
% At each step k the input u_k is now the TIME-AVERAGE of the flatness-
% consistent continuous input over the interval:
%
%       u_k = (1/Ts) * integral_{t_k}^{t_k+Ts} u(t) dt
%
% computed EXACTLY via an auxiliary ODE:
%
%       dz/dt = u(t),   z(t_k) = 0   =>   u_k = z(t_k+Ts) / Ts
%
% integrated with ode45.
%
% INPUTS
%   traj_obj     - trajectory object  (ShapeCircle | ShapeSpiral , ...)
%   T_sim        - total simulation duration  [s]
%   Ts           - sampling time              [s]
%   model        - nonlinear model object
%   flatness     - flatness map object (related to the nonlinear model)
%   options      - ode45 options struct (odeset) for the STATE propagation
%   quad_options - (optional) ode45 options struct (odeset) for the
%                  auxiliary u(t)-quadrature ODE. Default:
%                  odeset('RelTol',1e-10,'AbsTol',1e-12) — tight, since
%                  cost is irrelevant offline.
%
% OUTPUTS
%   x_ref     - nx x N_steps  state   trajectory matrix (nx = 12 or 13)
%   u_ref     -  4 x N_steps  input   trajectory matrix (now interval-averaged)
%   t_steps   -  1 x N_steps  time    vector  [s]

    %% Optional argument
    if nargin < 7 || isempty(quad_options)
        quad_options = odeset('RelTol', 1e-10, 'AbsTol', 1e-12);
    end

    %% Model parameters
    mp = ModelParameters();

    %% Time grid
    N_steps = round(T_sim / Ts) + 1;
    t_steps = linspace(0, T_sim, N_steps);   

    %% Initial condition via flatness map at t = 0 
    [s0, ds0, dds0, ddds0, dddds0] = traj_obj.get_flat_outputs(0);
    [x0, ~] = flatness.map(mp, s0, ds0, dds0, ddds0, dddds0);

    % State dimension (nx = 12 for Euler, 13 for Quaternion)
    nx = length(x0);
    if nx == 13
        % Normalize initial quaternion to prevent drift
        x0(7:10) = x0(7:10) / norm(x0(7:10));
    end

    %% Pre-allocate outputs
    x_ref = zeros(nx, N_steps);
    u_ref = zeros(4,  N_steps);

    x_ref(:, 1) = x0;
    x_curr      = x0;

    %% PWC Inputs integration
    for k = 1 : N_steps - 1

        t_now = t_steps(k);

        % Exact time-averaged input over [t_now, t_now+Ts] via auxiliary ODE
        u_k = computeMeanInput_ode45(traj_obj, mp, flatness, t_now, Ts, quad_options);

        % Store (now averaged, not point-sampled) input at step k
        u_ref(:, k) = u_k;

        % Integrate nonlinear dynamics with u_k held constant (unchanged)
        tspan = [t_now, t_now + Ts];
        [~, x_sol] = ode45(@(t, x) model.dynamics(t, x, u_k), ...
                           tspan, x_curr, options);

        x_curr = x_sol(end, :)';
        
        % Normalize quaternion for the 13-state representation
        if nx == 13
            x_curr(7:10) = x_curr(7:10) / norm(x_curr(7:10));
        end
        
        x_ref(:, k+1) = x_curr;

    end

    % Unwrap yaw (12 states)
    if nx == 12
        x_ref(9,:) = unwrap(x_ref(9,:));    
    end

    %% Repeat last input to fill the final column
    u_ref(:, N_steps) = u_ref(:, N_steps - 1);

end


function u_bar = computeMeanInput_ode45(traj_obj, mp, flatness, t_now, Ts, quad_options)
% COMPUTEMEANINPUT_ODE45  Exact (to solver tolerance) time-average of
%                         u(t) = flatness.map(...) over [t_now, t_now+Ts],
%                         via the auxiliary ODE  dz/dt = u(t), z(t_now)=0.
%
%   u_bar = z(t_now+Ts) / Ts

    z0    = zeros(4, 1);
    tspan = [t_now, t_now + Ts];

    [~, z_sol] = ode45(@(t, z) evalU(t, traj_obj, mp, flatness), ...
                        tspan, z0, quad_options);

    u_bar = z_sol(end, :)' / Ts;

end


function u = evalU(t, traj_obj, mp, flatness)
% EVALU  Evaluates the flatness-consistent input u(t) at a single time
%        instant. Used as the (state-independent) vector field of the
%        auxiliary quadrature ODE dz/dt = u(t).

    [s, ds, dds, ddds, dddds] = traj_obj.get_flat_outputs(t);
    [~, u] = flatness.map(mp, s, ds, dds, ddds, dddds);

end