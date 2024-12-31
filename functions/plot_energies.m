function [] = plot_energies(n, l, m, time_array, KE_radial_traj, KE_theta_traj, KE_phi_traj, V_traj, E_traj)
    hold on;
    
    theor_E_r = expected_E_r(n, l);
    if ~isnan(theor_E_r)
        yline(theor_E_r, 'LineStyle', '-', 'Color', [0, 0, 1], 'LineWidth', 0.7, 'HandleVisibility', 'off');
    end
    plot(time_array, KE_radial_traj, ...
        'LineStyle', '-', 'Color', [0, 0, 1], 'LineWidth', 2, 'DisplayName', '$\left< E_r \right>$');
    
    theot_E_theta = expected_E_theta(l, m);
    yline(theot_E_theta, 'LineStyle', '-', 'Color', [0, 1, 0], 'LineWidth', 0.7, 'HandleVisibility', 'off');
    plot(time_array, KE_theta_traj, ...
        'LineStyle', '-', 'Color', [0, 1, 0], 'LineWidth', 2, 'DisplayName', '$\left< E_\theta \right>$');

    if m ~= 0 % no need to plot trivial energy E_phi=0
        theot_E_phi = expected_E_phi(l, m);
        yline(theot_E_phi, 'LineStyle', '-', 'Color', [0, 1, 0.7], 'LineWidth', 0.7, 'HandleVisibility', 'off');
        plot(time_array, KE_phi_traj, ...
            'LineStyle', '-', 'Color', [0, 1, 0.7], 'LineWidth', 2, 'DisplayName', '$\left< E_\theta \right>$');
    end

    theor_V = expected_V(n);
    yline(theor_V, 'LineStyle', '-', 'Color', [1, 0.5, 0], 'LineWidth', 0.5, 'HandleVisibility', 'off');
    plot(time_array, V_traj, ...
        'LineStyle', '-', 'Color', [1, 0.5, 0], 'LineWidth', 2, 'DisplayName', '$\left< V \right>$');
   
    theor_E = expected_E(n);
    yline(theor_E, 'LineStyle', '-', 'Color', [1, 0, 0], 'LineWidth', 0.7, 'HandleVisibility', 'off');
    plot(time_array, E_traj, ...
        'LineStyle', '-', 'Color', [1, 0, 0], 'LineWidth', 2, 'DisplayName', '$\left< E \right>$');
    hold off;
    ylabel('Average Energies (J)');
    xlabel('Time (s)');
    legend('Interpreter', 'latex');
    
    if length(time_array) > 1
        xlim([0, time_array(end)]);
    end
    
    % Optimized for 2p0 and 2s0 states. Modify for other states.
    if (n == 2 && l == 1 && m == 0) || (n == 2 && l == 0 && m == 0)
        yticks(-15e-19:3e-19:15e-19);
        ylim([-15e-19, 15e-19]);
    end

end

function [exp_E_r] = expected_E_r(n, l) 
    [hbar, m_e, a_0, ~, ~, ~] = constants();
    mean_inv_r2 = 1/ (12 * a_0^2);

    exp_E_r = hbar^2 / ( 2 * m_e) * (1 / (a_0^2 * n^2) - mean_inv_r2 * l * (l + 1));
end

function [exp_E_theta] = expected_E_theta(l, m)
    [hbar, m_e, a_0, ~, ~, ~] = constants();
    mean_inv_r2 = 1/ (12 * a_0^2);

    exp_E_theta = hbar^2 / (2 * m_e) * mean_inv_r2 * (2 * l^2 - 2 * l * (m - 1) - m) / 2;
end

function [exp_E_phi] = expected_E_phi(l, m)
    [hbar, m_e, a_0, ~, ~, ~] = constants();
    mean_inv_r2 = 1/ (12 * a_0^2);

    exp_E_phi = hbar^2 / (2 * m_e) * mean_inv_r2 * (2 * l + 1) * m / 2;
end

function [exp_V] = expected_V(n)
    % Using the Virial Theorem:
    exp_V  = 2 * expected_E(n);
end

function [exp_E] = expected_E(n)
    % expected_total_energy = -5.4494e-19; l=1, m=0 
    [hbar ,m_e, a_0, ~, ~, ~] = constants();

    %expected_E  = - (m_e * e_charge^4) / ( 32 * pi^2 * epsilon_0^2 * hbar^2 * n^2);
    exp_E  = - hbar^2/( 2 * m_e) / (a_0^2 * n^2);
end