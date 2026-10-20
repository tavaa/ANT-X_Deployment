function [x_ref, u_ref, t_steps] = B01_GeneratePWC_reference(traj_obj, T_sim, Ts, options, quad_options)
% B01_GENERATEPWC_REFERENCE  Generates a Piece-Wise Constant (PWC)
%                             reference for the B01 model, by
%                             integrating the nonlinear NO-DELAY model
%                             forward with the interval-averaged input
%                             held constant at each sampling instant.
%
%     Flatness map : B01_FlatnessMap (always returns a 16-element state;
%     truncated to 12 here, since motor delay is not invertible via
%     flatness)
%     Model         : B01_ANTX_quadcopter(mp, false), i.e. delay_motors
%     FORCED to false.
%
% At each step k, u_k is the TIME-AVERAGE of the flatness-consistent
% continuous input over [k*Ts, (k+1)*Ts]:
%
%     u_k = (1/Ts) * integral_{t_k}^{t_k+Ts} u(t) dt
%
% computed exactly (to solver tolerance) via an auxiliary ODE
% dz/dt = u(t). 
%
% INPUTS
%   traj_obj     - trajectory object (ShapeCircle | ShapeSpiral | ...)
%   T_sim        - total simulation duration  [s]
%   Ts           - sampling time              [s]
%   options      - ode45 options struct (odeset) for the STATE propagation
%   quad_options - (optional) ode45 options struct (odeset) for the
%                  auxiliary u(t)-quadrature ODE. Default:
%                  odeset('RelTol',1e-10,'AbsTol',1e-12) -- tight, since
%                  cost is irrelevant offline.
%
% OUTPUTS
%   x_ref     - 12 x N_steps  state   trajectory matrix (rigid body only,
%               no motor-delay states -- see header note above)
%   u_ref     -  4 x N_steps  input   trajectory matrix [T1;T2;T3;T4]
%   t_steps   -  1 x N_steps  time    vector  [s]

    %% Optional argument
    if nargin < 5 || isempty(quad_options)
        quad_options = odeset('RelTol', 1e-10, 'AbsTol', 1e-12);
    end

    %% Model parameters and model instance (no-delay)
    mp    = ModelParameters();
    model = B01_ANTX_quadcopter(mp, false);   % delay_motors FORCED false

    nx = 12;   % fixed for this no-delay model

    %% Time grid
    N_steps = round(T_sim / Ts) + 1;
    t_steps = linspace(0, T_sim, N_steps);

    %% Initial condition via flatness map at t = 0 
    [s0, ds0, dds0, ddds0, dddds0] = traj_obj.get_flat_outputs(0);
    [x0_full, ~] = B01_FlatnessMap.map(mp, s0, ds0, dds0, ddds0, dddds0);
    x0 = x0_full(1:nx);

    %% Pre-allocate outputs
    x_ref = zeros(nx, N_steps);
    u_ref = zeros(4,  N_steps);

    x_ref(:, 1) = x0;
    x_curr      = x0;

    %% PWC input integration
    for k = 1 : N_steps - 1

        t_now = t_steps(k);

        % Exact time-averaged input over [t_now, t_now+Ts] via auxiliary ODE
        u_k = computeMeanInput_ode45(traj_obj, mp, t_now, Ts, quad_options);

        % Store (averaged) input at step k
        u_ref(:, k) = u_k;

        % Integrate no-delay nonlinear dynamics with u_k held constant
        tspan = [t_now, t_now + Ts];
        [~, x_sol] = ode45(@(t, x) model.dynamics(t, x, u_k), ...
                           tspan, x_curr, options);

        x_curr = x_sol(end, :)';
        x_ref(:, k+1) = x_curr;

    end

    % Unwrap yaw 
    x_ref(9,:) = unwrap(x_ref(9,:));

    %% Repeat last input to fill the final column
    u_ref(:, N_steps) = u_ref(:, N_steps - 1);

end


function u_bar = computeMeanInput_ode45(traj_obj, mp, t_now, Ts, quad_options)
% COMPUTEMEANINPUT_ODE45  Exact (to solver tolerance) time-average of
%                         u(t) = B01_FlatnessMap.map(...) over
%                         [t_now, t_now+Ts], via the auxiliary ODE
%                         dz/dt = u(t), z(t_now)=0.
%
%   u_bar = z(t_now+Ts) / Ts

    z0    = zeros(4, 1);
    tspan = [t_now, t_now + Ts];

    [~, z_sol] = ode45(@(t, z) evalU(t, traj_obj, mp), ...
                        tspan, z0, quad_options);

    u_bar = z_sol(end, :)' / Ts;

end


function u = evalU(t, traj_obj, mp)
% EVALU  Evaluates the flatness-consistent input u(t) = [T1;T2;T3;T4] at
%        a single time instant. Used as the (state-independent) vector
%        field of the auxiliary quadrature ODE dz/dt = u(t).
%
%        Only uref is needed here (not xref) 

    [s, ds, dds, ddds, dddds] = traj_obj.get_flat_outputs(t);
    [~, u] = B01_FlatnessMap.map(mp, s, ds, dds, ddds, dddds);

end