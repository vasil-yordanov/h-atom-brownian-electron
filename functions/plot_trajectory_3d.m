function [] = plot_trajectory_3d(n, l, r_traj, theta_traj, phi_traj, a_0, traj_points, params_set_name, show_all, trajectory_style)
    % Set default for show_all based on the number of input arguments
    if nargin < 9 || isempty(show_all)
        show_all = false;
    end
    if nargin < 10 || isempty(trajectory_style)
        trajectory_style = 'radius';
    end

    % Determine the number of trajectory points to process
    m = min(traj_points, size(r_traj, 2));
    if m < 2
        return;
    end

    % Determine indices based on show_all flag
    if ~show_all
        % Fixed decimation stride based on the FULL trajectory length (constant),
        % not the growing traj_points. A traj_points-dependent stride changes the
        % sampled points every frame and makes the trace jump/swim; a fixed
        % stride keeps the vertices stable so the trace just grows by appending.
        every = max(floor(size(r_traj, 2) / 1e6), 1);
        idx = 1:every:m;
    else
        idx = 1:m;
    end
    
    % Precompute repeated terms to avoid redundant calculations
    r = r_traj(idx);
    theta = theta_traj(idx);
    phi = phi_traj(idx);
    
    s_theta = sin(theta);
    c_phi = cos(phi);
    s_phi = sin(phi);
    c_theta = cos(theta);
    
    % Compute Cartesian coordinates in a_0 units, and colour the line by radius.
    x = r .* s_theta .* c_phi ./ a_0;
    y = r .* s_theta .* s_phi ./ a_0;
    z = r .* c_theta ./ a_0;
    c = r ./ a_0;

    ax = gca;
    if strcmp(trajectory_style, 'sphere')
        plot3(ax, x, y, z, 'Color', [0 0.4470 0.7410], ...
              'LineWidth', 0.8);
        colorbar(ax, 'off');
    else
        surface(ax, [x; x], [y; y], [z; z], [c; c], ...
                'FaceColor', 'none', 'EdgeColor', 'interp', ...
                'LineWidth', 0.7, 'EdgeAlpha', 0.7);
        colormap(ax, turbo);
        colorbar(ax, 'off');
    end
    xlabel(ax, 'X / a_0', 'FontSize', 10);
    ylabel(ax, 'Y / a_0', 'FontSize', 10);
    zlabel(ax, 'Z / a_0', 'FontSize', 10);

    % Cubic axis box, shared with export_cloud_full via cloud_cube_limit so the
    % manuscript clouds and the live/movie clouds use the same geometry. The
    % 1s0_1ps cube is snug (its drift tail is intentionally cropped); all other
    % states use half the histogram range, padded. Set the limits FIRST.
    lim = cloud_cube_limit(n, l, a_0, params_set_name);
    xlim([-lim lim]);
    ylim([-lim lim]);
    zlim([-lim lim]);

    % Pin the box geometry to FIXED values so it is identical every frame. The
    % limits above are a cube, so DataAspectRatio = PlotBoxAspectRatio = [1 1 1]
    % renders an undistorted cube; making both manual also disables MATLAB's
    % stretch-to-fill (which would shrink the rotated cube and leave a big white
    % margin). Do NOT use `axis equal; axis vis3d;` here: vis3d freezes the
    % plot-box aspect ratio from the *current data extent*, which grows every
    % frame, so the box would rescale / jump from frame to frame in the movie.
    clim(ax, [0 lim]);
    set(ax, 'DataAspectRatio', [1 1 1], 'PlotBoxAspectRatio', [1 1 1], 'FontSize', 9);
    camup([0 0 1]);          % fix +z as "up"
    if strcmp(params_set_name, '1s0_1ps')
        view(-18, 16);
    else
        view(-15, 10);
    end

    grid on;

end
