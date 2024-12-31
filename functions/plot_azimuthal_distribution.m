function [] = plot_azimuthal_distribution(hist_bins_phi, hist_counts_phi, phi_values, P_phi_analytic)
    % Normalize the Histogram Counts
    hist_counts_phi_norm = normalize_histogram(hist_counts_phi, hist_bins_phi);

    % Calculate bin centers for phi histogram
    phi_bin_centers = hist_bins_phi(1:end-1) + diff(hist_bins_phi)/2;

    bar(phi_bin_centers, hist_counts_phi_norm, 'FaceAlpha', 0.6, 'EdgeColor', 'none');
    hold on;
    plot(phi_values, P_phi_analytic, 'r', 'LineWidth', 2);
    hold off;
    xlabel('$\phi$ (radians)', 'Interpreter', 'latex');
    ylabel('Azimuthal $(\phi)$ Distribution', 'Interpreter', 'latex');
    legend('Simulation', 'Analytic Solution');
    ylim([0, max(max(P_phi_analytic), max(hist_counts_phi_norm))*1.25]);
end

