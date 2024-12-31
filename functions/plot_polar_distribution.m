function [] = plot_polar_distribution(hist_bins_theta, hist_counts_theta, theta_values, P_theta_analytic)
    % Normalize the Histogram Counts
    hist_counts_theta_norm = normalize_histogram(hist_counts_theta, hist_bins_theta);

    % Calculate bin centers for theta histogram
    theta_bin_centers = hist_bins_theta(1:end-1) + diff(hist_bins_theta)/2;
 
    bar(theta_bin_centers, hist_counts_theta_norm, 'FaceAlpha', 0.6, 'EdgeColor', 'none');
    hold on;
    plot(theta_values, P_theta_analytic, 'r', 'LineWidth', 2);
    hold off;
   
    xlabel('$\theta$ (radians)', 'Interpreter', 'latex');
    ylabel('Polar $(\theta)$ Distribution', 'Interpreter', 'latex');
    
    legend('Simulation', 'Analytical Solution');

    ylim([0, max(max(P_theta_analytic), max(hist_counts_theta_norm))*1.25]);
end

