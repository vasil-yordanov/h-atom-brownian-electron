function [] = plot_radial_distribution(hist_bins_r, hist_counts_r, r_values, P_analytic_r)
    % Bohr radius
    a_0 = 5.29177e-11;

    % Normalize the Histogram Counts
    hist_counts_r_norm = normalize_histogram(hist_counts_r, hist_bins_r);
 
    % Calculate bin centers for R histogram
    r_bin_centers = hist_bins_r(1:end-1) + diff(hist_bins_r)/2;       

    bar(r_bin_centers / a_0, hist_counts_r_norm, 'FaceAlpha', 0.6, 'EdgeColor', 'none');
    hold on;
    plot(r_values / a_0, P_analytic_r, 'r', 'LineWidth', 2);
    hold off;
    xlabel('$r / a_0$', 'Interpreter', 'latex');
    ylabel('Radial $(r)$ Distribution', 'Interpreter', 'latex');
    legend('Simulation', 'Analytical Solution');
    ylim([0, max(max(P_analytic_r), max(hist_counts_r_norm))*1.25]);
end

