function [] = plot_trajectory_xy(n, l, r_traj, theta_traj, phi_traj, a_0, traj_points, params_set_name, show_all)
    if nargin == 8
        show_all = false;
    end

    m = min(traj_points, size(r_traj, 2));
    if m < 2
        return;
    end

    if ~show_all
        every = max(floor(traj_points / 1e5), 1);
        idx = 1:every:m;
    else
        idx = 1:m;
    end

    r = r_traj(idx);
    theta = theta_traj(idx);
    phi = phi_traj(idx);

    x = r .* sin(theta) .* cos(phi) ./ a_0;
    y = r .* sin(theta) .* sin(phi) ./ a_0;

    c = linspace(0, 1, numel(x));
    surface([x; x], [y; y], zeros(2, numel(x)), [c; c], ...
        'FaceColor', 'none', 'EdgeColor', 'interp', 'LineWidth', 1.6);
    blue_gamma = [linspace(0.82, 0.05, 256)', linspace(0.90, 0.20, 256)', linspace(1.00, 0.78, 256)'];
    colormap(gca, blue_gamma);
    cb = colorbar('eastoutside');
    cb.Label.String = 'Normalized time';
    cb.FontSize = 12;
    cb.Label.FontSize = 14;
    caxis([0, 1]);

    xlabel('X', 'FontSize', 14);
    ylabel('Y', 'FontSize', 14);
    set(gca, 'FontSize', 12);
    axis equal;
    grid on;

    [~, r_max] = radial_histogram_range(n, l, a_0, params_set_name);
    if strcmp(params_set_name, '1s0_1ps')
        r_max = 5 * a_0;
    else
        r_max = 0.5 * r_max;
    end
    xlim([-1.5 * r_max / a_0, 1.5 * r_max / a_0]);
    ylim([-1.5 * r_max / a_0, 1.5 * r_max / a_0]);
end
