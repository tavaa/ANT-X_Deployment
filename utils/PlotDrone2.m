classdef PlotDrone2
% Coordinate conventions:
%   Inertial frame : NED  (x North, y East, z Down)
%   Body frame     : FRD  (x Forward, y Right, z Down)
%   Euler sequence : Z-Y-X (yaw-pitch-roll)
%
% Sign conventions (consistent with FRD + NED):
%   phi   > 0  -> right wing DOWN     (positive roll)
%   theta > 0  -> nose UP           (positive pitch)
%   psi   > 0  -> clockwise yaw viewed from above
%
% Plot convention:
%   ZDir, YDir = reverse is used here.
%   Therefore:
%       z_NED = 0      -> ground plane
%       z_NED < 0      -> altitude above ground
%       z_NED > 0      -> below ground
%
% Drone geometry:
%   X-configuration quadrotor.
%
%   Body axes relative to motors:
%
%                      x_B (Forward)
%                             ^
%                             |
%                      M1     |     M4
%                         \   |   /
%                          \  |  /
%           left  <--------- CoG ---------> right
%                          /  |  \
%                         /   |   \
%                      M2     |     M3
%
%                             v
%
%   x_B passes BETWEEN:
%       front motors : M1-M4
%       rear motors  : M2-M3
%
%   y_B passes BETWEEN:
%       right motors : M4-M3
%       left motors  : M1-M2
%
% Motor numbering used in THIS code (Ghignoni convention):
%   M1 : front-left   (+x_B,-y_B)
%   M2 : rear-left    (-x_B,-y_B)
%   M3 : rear-right   (-x_B,+y_B)
%   M4 : front-right  (+x_B,+y_B)
%
% Usage:
%   h = PlotDrone2.init(ax, mp);
%   PlotDrone2.update(h, pos_NED, R);

    properties (Constant)

        % Body axis colors
        COL_XB = [0.85, 0.15, 0.15];   % red   (x_B Forward)
        COL_YB = [0.20, 0.70, 0.20];   % green (y_B Right)
        COL_ZB = [0.15, 0.35, 0.85];   % blue  (z_B Down)
        COL_NED_X = [0.9, 0.0, 0.0];
        COL_NED_Y = [0.0, 0.6, 0.0];
        COL_NED_Z = [0.0, 0.0, 0.9];

        % Drone geometry [m]
        FRAME_HX = 0.10;   % half-length x body (total 20cm)
        FRAME_HY = 0.10;   % half-width  y body (total 20cm)
        FRAME_HZ = 0.02;   % half-height z body (total 4cm)

        % Visualization scales
        ARROW_LEN  = 0.10;   % body axis arrow length [m]
        NED_ARROW  = 0.30;   % NED reference arrow length [m]
        MOTOR_R    = 0.018;  % motor disk radius [m]
        COG_R      = 6;      % CoG marker size [pt]
        MOTOR_SZ   = 8;      % motor circle size [pt]
        TRAJ_SKIP  = 5;      % animate every N steps

    end

    methods (Static)

        function setup_axes(ax, title_str)
        %% SETUP_AXES  Configure NED plot axes with correct z orientation.
        %
        %   With ZDir = 'reverse': more negative z_NED (higher altitude)
        %   appears visually higher in the plot. Ground at z=0 is at
        %   the visual bottom.

            axes(ax);
            set(ax, 'ZDir', 'reverse');
            set(ax, 'YDir', 'reverse');
            %set(ax, 'XDir', 'reverse');
            grid(ax, 'on');
            axis(ax, 'equal');
            xlabel(ax, 'x_I  North [m]');
            ylabel(ax, 'y_I  East  [m]');
            zlabel(ax, 'z_{NED}  (neg = up) [m]');
            if nargin > 1
                title(ax, title_str);
            end
            view(ax, [-35, 25]);
        end


        function draw_ground(ax, x_range, y_range, z_ground)
        %% DRAW_GROUND  Draw ground plane at z_ground (default 0 in NED).
            if nargin < 4, z_ground = 0; end
            [xg, yg] = meshgrid(x_range, y_range);
            zg = z_ground * ones(size(xg));
            surf(ax, xg, yg, zg, 'FaceAlpha',0.15, 'FaceColor',[0.6 0.8 0.4], ...
                'EdgeColor',[0.5 0.7 0.3], 'LineWidth',0.5);
        end

        function h = init(ax, mp)
        %% INIT  Initialize all graphic handles for the drone.
        %
        %   ax : axes handle
        %   mp : ModelParameters (needs fields b, m)
        %
        %   Returns struct h with all graphic handles.

            hold(ax, 'on');

            b = mp.b / sqrt(2);   % arm projection onto diagonal [m]

            % motor positions in body frame (X-config, 45 deg)
            %   M1 : front-left   (+x_B,-y_B)
            %   M2 : rear-left    (-x_B,-y_B)
            %   M3 : rear-right   (-x_B,+y_B)
            %   M4 : front-right  (+x_B,+y_B)
            motors_b = [ ...
                b, -b, -b,  b;
                -b, -b,  b,  b;
                0,  0,  0,  0];

            % frame corners in body frame (8 corners of box)
            fx = PlotDrone2.FRAME_HX;
            fy = PlotDrone2.FRAME_HY;
            fz = PlotDrone2.FRAME_HZ;
            corners_b = [ fx,  fx,  fx,  fx, -fx, -fx, -fx, -fx;
                           fy,  fy, -fy, -fy,  fy,  fy, -fy, -fy;
                           fz, -fz,  fz, -fz,  fz, -fz,  fz, -fz];

            % 12 edges of the box (pairs of corner indices)
            edges = [1,2; 3,4; 5,6; 7,8;   % longitudinal edges
                     1,3; 2,4; 5,7; 6,8;   % lateral edges
                     1,5; 2,6; 3,7; 4,8];  % vertical edges

            % initialize arm lines (4 arms: CoG to each motor) 
            h.arms = gobjects(4, 1);
            arm_colors = {'k','k','k','k'};
            for i = 1:4
                h.arms(i) = plot3(ax, NaN,NaN,NaN, '-', ...
                    'Color', arm_colors{i}, 'LineWidth', 2.5);
            end

            % initialize frame box edges (dashed) 
            h.frame_edges = gobjects(12, 1);
            for e = 1:12
                h.frame_edges(e) = plot3(ax, NaN,NaN,NaN, 'k--', ...
                    'LineWidth', 0.8);
            end

            % CoG (filled black circle) 
            h.cog = plot3(ax, NaN, NaN, NaN, 'ko', ...
                'MarkerSize', PlotDrone2.COG_R, 'MarkerFaceColor', 'k');

            % Motors (hollow circles, drawn as small rings) 
            h.motors = gobjects(4, 1);
            for i = 1:4
                h.motors(i) = plot3(ax, NaN,NaN,NaN, 'o', ...
                    'MarkerSize', PlotDrone2.MOTOR_SZ, ...
                    'MarkerEdgeColor', [0.2 0.2 0.2], ...
                    'MarkerFaceColor', 'w', 'LineWidth', 1.5);
            end

            % Motor numbers
            h.motor_labels = gobjects(4,1);

            for i = 1:4
                h.motor_labels(i) = text(ax, NaN, NaN, NaN, sprintf('%d',i), ...
                'HorizontalAlignment','center', ...
                'VerticalAlignment','middle', ...
                'FontWeight','bold', ...
                'FontSize',8, ...
                'Color',[0 0 0]);
            end

            % Motor disks (thin circles in body xy plane) 
            h.disks = gobjects(4, 1);
            for i = 1:4
                h.disks(i) = plot3(ax, NaN,NaN,NaN, 'Color', [0.5 0.5 0.5], ...
                    'LineWidth', 0.8);
            end

            % FRD body axis arrows 
            L = PlotDrone2.ARROW_LEN;
            h.xb_arr = quiver3(ax, 0,0,0, L,0,0, 0, ...
                'Color', PlotDrone2.COL_XB, 'LineWidth', 2.0, 'MaxHeadSize', 0.5);
            h.yb_arr = quiver3(ax, 0,0,0, 0,L,0, 0, ...
                'Color', PlotDrone2.COL_YB, 'LineWidth', 2.0, 'MaxHeadSize', 0.5);
            h.zb_arr = quiver3(ax, 0,0,0, 0,0,L, 0, ...
                'Color', PlotDrone2.COL_ZB, 'LineWidth', 2.0, 'MaxHeadSize', 0.5);

            % FRD axis labels 
            h.xb_lbl = text(ax, NaN,NaN,NaN, 'x_B (F)', ...
                'Color', PlotDrone2.COL_XB, 'FontSize', 7, 'FontWeight','bold');
            h.yb_lbl = text(ax, NaN,NaN,NaN, 'y_B (R)', ...
                'Color', PlotDrone2.COL_YB, 'FontSize', 7, 'FontWeight','bold');
            h.zb_lbl = text(ax, NaN,NaN,NaN, 'z_B (D)', ...
                'Color', PlotDrone2.COL_ZB, 'FontSize', 7, 'FontWeight','bold');

            % legend patch (invisible anchor for legend) 
            h.leg_xb = plot3(ax,NaN,NaN,NaN,'-','Color',PlotDrone2.COL_XB,'LineWidth',2);
            h.leg_yb = plot3(ax,NaN,NaN,NaN,'-','Color',PlotDrone2.COL_YB,'LineWidth',2);
            h.leg_zb = plot3(ax,NaN,NaN,NaN,'-','Color',PlotDrone2.COL_ZB,'LineWidth',2);
            h.leg_arm = plot3(ax,NaN,NaN,NaN,'k-','LineWidth',2.5);

            legend(ax, [h.leg_xb, h.leg_yb, h.leg_zb, h.leg_arm], ...
                {'x_B Forward (FRD)','y_B Right (FRD)','z_B Down (FRD)','Arms (X-config)'}, ...
                'Location','northeast','FontSize',7);

            % store body geometry for use in update
            h.motors_b = motors_b;
            h.corners_b = corners_b;
            h.edges = edges;
            h.b = b;

        end

        function update(h, pos_NED, R)
        %% UPDATE  Redraw drone at given NED position and rotation.
        %
        %   pos_NED : [3x1] position in NED [m]
        %   R       : [3x3] rotation matrix body->NED (DCM)
        %
        % The drone body axes x_B, y_B, z_B are the columns of R.
        %   R(:,1) = x_B in NED  (Forward)
        %   R(:,2) = y_B in NED  (Right)
        %   R(:,3) = z_B in NED  (Down)

            p = pos_NED(:);   % ensure column

            % Rotate geometry into NED
            motors_NED  = R * h.motors_b + p;
            corners_NED = R * h.corners_b + p;

            % Arms (CoG to each motor)
            for i = 1:4
                set(h.arms(i), ...
                    'XData', [p(1), motors_NED(1,i)], ...
                    'YData', [p(2), motors_NED(2,i)], ...
                    'ZData', [p(3), motors_NED(3,i)]);
            end

            % Frame box edges
            for e = 1:12
                ci = h.edges(e,1);
                cj = h.edges(e,2);
                set(h.frame_edges(e), ...
                    'XData', [corners_NED(1,ci), corners_NED(1,cj)], ...
                    'YData', [corners_NED(2,ci), corners_NED(2,cj)], ...
                    'ZData', [corners_NED(3,ci), corners_NED(3,cj)]);
            end

            % CoG
            set(h.cog, 'XData',p(1), 'YData',p(2), 'ZData',p(3));

            % Motors (hollow circles at arm tips)
            for i = 1:4
                set(h.motors(i), ...
                    'XData', motors_NED(1,i), ...
                    'YData', motors_NED(2,i), ...
                    'ZData', motors_NED(3,i));
            end

            % Motor numbers
            for i = 1:4
                set(h.motor_labels(i), ...
                'Position', motors_NED(:,i));
            end

            % Motor disks (small circles in body xy plane)
            n_pts = 24;
            ang = linspace(0, 2*pi, n_pts);
            r_disk = PlotDrone2.MOTOR_R;
            for i = 1:4
                % circle in body xy plane centered at motor i
                circle_b = h.motors_b(:,i) + ...
                    [r_disk*cos(ang); r_disk*sin(ang); zeros(1,n_pts)];
                circle_NED = R * circle_b + p;
                set(h.disks(i), ...
                    'XData', circle_NED(1,:), ...
                    'YData', circle_NED(2,:), ...
                    'ZData', circle_NED(3,:));
            end

            % FRD body axis arrows
            L = PlotDrone2.ARROW_LEN;
            xb = R(:,1);  yb = R(:,2);  zb = R(:,3);

            set(h.xb_arr, 'XData',p(1),'YData',p(2),'ZData',p(3), ...
                'UData',L*xb(1),'VData',L*xb(2),'WData',L*xb(3));
            set(h.yb_arr, 'XData',p(1),'YData',p(2),'ZData',p(3), ...
                'UData',L*yb(1),'VData',L*yb(2),'WData',L*yb(3));
            set(h.zb_arr, 'XData',p(1),'YData',p(2),'ZData',p(3), ...
                'UData',L*zb(1),'VData',L*zb(2),'WData',L*zb(3));

            % FRD axis labels (at arrow tips)
            off = 1.15 * L;
            set(h.xb_lbl, 'Position', p + off*xb);
            set(h.yb_lbl, 'Position', p + off*yb);
            set(h.zb_lbl, 'Position', p + off*zb);

        end

    end
end