function lim = cloud_cube_limit(n, l, a_0, params_set_name)
    % CLOUD_CUBE_LIMIT  Half-width (in units of a_0) of the cubic axis box for a 3D
    % trajectory cloud. Shared by plot_trajectory_3d (live plot / dashboard movie)
    % and export_cloud_full (manuscript figures) so the box is defined in one place.
    % The box is centred at the origin: the axis limits are [-lim, lim] on x, y, z.
    [~, r_max] = radial_histogram_range(n, l, a_0, params_set_name);
    if strcmp(params_set_name, '1s0_1ps')
        % Snug cube cropping the 1s drift tail (the electron starts at r = 10 a_0);
        % a wider box would leave large vertical padding.
        lim = 1.2 * (5 * a_0) / a_0;
    else
        lim = 1.5 * (r_max * 0.5) / a_0;
    end
end
