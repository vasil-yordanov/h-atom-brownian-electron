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
        % Calculate the step size, ensuring it's at least 1
        every = max(floor(traj_points / 1e6), 1);
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

    plot3(x ./ a_0, y ./ a_0, z ./ a_0, '-');
    xlabel('X');
    ylabel('Y');
    zlabel('Z');
    
    [~, r_max] = radial_histogram_range(n, l, a_0, params_set_name);
    
    if strcmp(params_set_name, '1s0_1fs')
        view(-15, 10);
    elseif strcmp(params_set_name, '1s0_1ps')
        r_max = 5 * a_0;
        xlim([-r_max/a_0 2*r_max/a_0]);
        ylim([-r_max/a_0 r_max/a_0]);
        zlim([-r_max/a_0 r_max/a_0]);
        yticks([-4 0 4]);
        view(-3, 15);
    else
        r_max = r_max * 0.5;
        xlim([-1.5*r_max/a_0 1.5*r_max/a_0]);
        ylim([-1.5*r_max/a_0 1.5*r_max/a_0]);
        zlim([-1.5*r_max/a_0 1.5*r_max/a_0]);
        view(-15, 10);
    end

    grid on;

end

