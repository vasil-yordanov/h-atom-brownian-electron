function [] = plot_azimuthal_current(time_array, b_phi_avg_array, n, l, m)
    [hbar, m_e, a_0, ~, ~, ~] = constants();

    hold on;
    yline(0, 'k--', 'LineWidth', 0.8, 'HandleVisibility', 'off');

    theory_handle = [];
    if m ~= 0
        theta_values = linspace(1e-6, pi - 1e-6, 20000);
        Y_theta = angular_wavefunction_theta(theta_values, l, m);
        angular_factor = 2 * pi * trapz(theta_values, abs(Y_theta).^2);
        radial_factor = 1 / (n^2 * a_0);
        b_phi_theory = (hbar / m_e) * m * radial_factor * angular_factor;
        theory_handle = yline(b_phi_theory, 'r--', 'LineWidth', 1.4);
    end

    mean_handle = plot(time_array, b_phi_avg_array, 'LineWidth', 1.8, 'Color', [0.2, 0.5, 0.9]);
    hold off;

    xlabel('Time (s)', 'FontSize', 14);
    ylabel('Mean Azimuthal Current - $\langle b_\phi \rangle$ (m/s)', 'Interpreter', 'latex', 'FontSize', 14);
    set(gca, 'FontSize', 12);

    if ~isempty(theory_handle)
        legend([mean_handle, theory_handle], {'Simulation', 'Theory'}, 'Location', 'best', 'FontSize', 11);
    end

    if length(time_array) > 1
        xlim([0, time_array(end)]);
    end
end
