function [] = plot_trajectory_3d(n, l, r_traj, theta_traj, phi_traj, a_0, traj_points, params_set_name, show_all)
    % Set default for show_all based on the number of input arguments
    if nargin == 8
        show_all = false;
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
    
    % Compute Cartesian coordinates
    x = r .* s_theta .* c_phi;
    y = r .* s_theta .* s_phi;
    z = r .* c_theta;

    plot3(x ./ a_0, y ./ a_0, z ./ a_0, '-', 'LineWidth', 1.0, 'Color', [0.1, 0.35, 0.8]);
    xlabel('X', 'FontSize', 14);
    ylabel('Y', 'FontSize', 14);
    zlabel('Z', 'FontSize', 14);

    % Rendering scale of the 3D cloud. Do NOT drop axis vis3d: without it MATLAB
    % stretch-to-fits the rotated plot box to the axes at every view angle, which
    % shrinks the cube and leaves a big empty margin -- exportgraphics then writes
    % a tiny cloud in a large white box. axis equal keeps the cube undistorted and
    % camup fixes +z as "up".
    axis equal;
    axis vis3d;
    camup([0 0 1]);
    set(gca, 'FontSize', 12);

    % Cubic axis box, shared with export_cloud_full via cloud_cube_limit so the
    % manuscript clouds and the live/movie clouds use the same geometry. The
    % 1s0_1ps cube is snug (its drift tail is intentionally cropped); all other
    % states use half the histogram range, padded.
    lim = cloud_cube_limit(n, l, a_0, params_set_name);
    xlim([-lim lim]);
    ylim([-lim lim]);
    zlim([-lim lim]);
    if strcmp(params_set_name, '1s0_1ps')
        view(-18, 16);
    else
        view(-15, 10);
    end

    grid on;

end

