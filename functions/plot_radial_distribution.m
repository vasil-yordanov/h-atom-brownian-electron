function [] = plot_radial_distribution(hist_bins_r, hist_counts_r, r_values, P_analytic_r)
    % Bohr radius
    a_0 = 5.29177e-11;

    % Normalize the histogram as a probability density in the plotted variable
    % r/a_0 (so the area under it is 1 and the values are O(1)), consistent with
    % the polar and azimuthal panels. Normalizing in r (metres) instead made the
    % density ~10^10 /m, so the y-axis read in units of 10^9.
    hist_counts_r_norm = normalize_histogram(hist_counts_r, hist_bins_r / a_0);

    % Calculate bin centers for R histogram
    r_bin_centers = hist_bins_r(1:end-1) + diff(hist_bins_r)/2;

    % P_analytic_r = r^2 |R_nl(r)|^2 is a density per metre (its integral over r
    % is 1); multiply by a_0 to express it per r/a_0, matching the histogram.
    P_analytic_r_norm = P_analytic_r * a_0;

    bar(r_bin_centers / a_0, hist_counts_r_norm, 'FaceAlpha', 0.6, 'EdgeColor', 'none');
    hold on;
    plot(r_values / a_0, P_analytic_r_norm, 'r', 'LineWidth', 2);
    hold off;
    xlabel('$r / a_0$', 'Interpreter', 'latex', 'FontSize', 14);
    ylabel('Radial $(r)$ Distribution', 'Interpreter', 'latex', 'FontSize', 14);
    legend('Simulation', 'Analytical Solution', 'FontSize', 11);
    set(gca, 'FontSize', 12);
    ylim([0, max(max(P_analytic_r_norm), max(hist_counts_r_norm))*1.25]);
end

