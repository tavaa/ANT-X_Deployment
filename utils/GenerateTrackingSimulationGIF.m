% Generates a GIF for MPC tracking

function GenerateTrackingSimulationGIF(log_states, x_ref_pwc, gif_filename, title_str, mode, mp, PlotDrone)

    % States
    X_sim = log_states(:, 2);
    Y_sim = log_states(:, 3);
    Z_sim = log_states(:, 4);

    phi_sim   = log_states(:, 8);
    theta_sim = log_states(:, 9);
    psi_sim   = log_states(:,10);

    N = length(X_sim);

    % References
    idx_end = min(N, size(x_ref_pwc, 2));

    X_ref = x_ref_pwc(1,1:idx_end);
    Y_ref = x_ref_pwc(2,1:idx_end);
    Z_ref = x_ref_pwc(3,1:idx_end);

    % Figure
    fig = figure('Color', 'w', ...
                 'Position', [100 100 1400 900]);

    ax = axes(fig);
    hold(ax, 'on');

    PlotDrone.setup_axes(ax, title_str);

    % Ground plane
    x_margin = 0.5;
    y_margin = 0.5;

    x_range = linspace(min([X_sim; X_ref']) - x_margin, ...
                       max([X_sim; X_ref']) + x_margin, 10);

    y_range = linspace(min([Y_sim; Y_ref']) - y_margin, ...
                       max([Y_sim; Y_ref']) + y_margin, 10);

    PlotDrone.draw_ground(ax, x_range, y_range, 0);

    if strcmp(mode, '2D')
        view(ax, 2);
    else
        view(ax, -35, 25);
    end

    % Reference Trajectory
    plot3(ax, ...
        X_ref, Y_ref, Z_ref, ...
        'k--', ...
        'LineWidth', 1.5, ...
        'DisplayName', 'Reference');

    % Simulated Trajectory
    h_traj = plot3(ax, ...
        NaN, NaN, NaN, ...
        'Color', [0 0.4 1], ...
        'LineWidth', 2.5, ...
        'DisplayName', 'Simulated');

    % Target point
    h_target = plot3(ax, ...
        NaN, NaN, NaN, ...
        'x', ...
        'MarkerSize', 14, ...
        'LineWidth', 2.5, ...
        'Color', [0.8 0.4 0]);

    % Drone Initializaion
    h_drone = PlotDrone.init(ax, mp);

    legend(ax, 'Location', 'best');

    axis(ax, 'padded');

    % GIF SETTINGS
    delay = 0.05;
    skipFrames = 1; % modify GIF duration here

    % ANIMATION LOOP
    for k = 1:skipFrames:N

        % Trajectory Update
        set(h_traj, ...
            'XData', X_sim(1:k), ...
            'YData', Y_sim(1:k), ...
            'ZData', Z_sim(1:k));

        idx = min(k, idx_end);

        set(h_target, ...
            'XData', X_ref(idx), ...
            'YData', Y_ref(idx), ...
            'ZData', Z_ref(idx));

        % CURRENT STATE
        pos = [
            X_sim(k);
            Y_sim(k);
            Z_sim(k)
        ];

        phi   = phi_sim(k);
        theta = theta_sim(k);
        psi   = psi_sim(k);


        % ROTATION MATRIX BODY -> NED
        % ZYX convention
        R = Euler2R_ZYX(phi, theta, psi);

        % UPDATE DRONE
        PlotDrone.update(h_drone, pos, R);
        drawnow;

        % FRAME EXPORT
        frame = getframe(fig);
        im = frame2im(frame);
        [A, map] = rgb2ind(im, 256);

        %% SAVE GIF
        %if k == 1

        %    imwrite(A, map, gif_filename, ...
        %        'gif', ...
        %        'LoopCount', inf, ...
        %        'DelayTime', delay);

        %else

        %    imwrite(A, map, gif_filename, ...
        %        'gif', ...
        %        'WriteMode', 'append', ...
        %        'DelayTime', delay);
        %end
    end

    close(fig);

end