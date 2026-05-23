function [] = plot_trajectory_3d(n, l, r_traj, theta_traj, phi_traj, a_0, traj_points, params_set_name, show_all)
    if nargin == 8
        show_all = false;
    end

    m = min(traj_points, size(r_traj, 2));
    if m < 2
        return;
    end

    if ~show_all
        every = max(floor(traj_points / 1e6), 1);
        idx = 1:every:m;
    else
        idx = 1:m;
    end

    r = r_traj(idx);
    theta = theta_traj(idx);
    phi = phi_traj(idx);

    x = r .* sin(theta) .* cos(phi);
    y = r .* sin(theta) .* sin(phi);
    z = r .* cos(theta);

    plot3(x ./ a_0, y ./ a_0, z ./ a_0, '-', 'LineWidth', 1.0, 'Color', [0.1, 0.35, 0.8]);
    xlabel('X', 'FontSize', 14);
    ylabel('Y', 'FontSize', 14);
    zlabel('Z', 'FontSize', 14);
    axis equal;
    axis vis3d;
    camup([0 0 1]);
    set(gca, 'FontSize', 12);

    [~, r_max] = radial_histogram_range(n, l, a_0, params_set_name);
    if strcmp(params_set_name, '1s0_1ps')
        r_max = 5 * a_0;
    else
        r_max = 0.5 * r_max;
    end

    xlim([-1.5 * r_max / a_0, 1.5 * r_max / a_0]);
    ylim([-1.5 * r_max / a_0, 1.5 * r_max / a_0]);
    zlim([-1.5 * r_max / a_0, 1.5 * r_max / a_0]);
    view(-18, 16);
    grid on;
end
