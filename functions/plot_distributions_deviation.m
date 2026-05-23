function [] = plot_distributions_deviation(time_array, deviation_array_R, deviation_array_theta, deviation_array_phi)
    hold on;
    plot(time_array, deviation_array_R, 'LineStyle', '-', 'Color', [0, 0, 1], 'LineWidth', 1.5);
    plot(time_array, deviation_array_theta, 'LineStyle', '-', 'Color', [1, 0, 0], 'LineWidth', 1.5);
    plot(time_array, deviation_array_phi, 'LineStyle', '-', 'Color', [0, 1, 0], 'LineWidth', 1.5);
    ylabel('Distribution Deviation (%)', 'FontSize', 14);
    hold off;

    xlabel('Time (s)', 'FontSize', 14);
    legend('Radial $(r)$', 'Polar $(\theta)$', 'Azimuthal $(\phi)$', 'Interpreter', 'latex', 'FontSize', 11);
    set(gca, 'FontSize', 12);

    if length(time_array) > 1
        xlim([0, time_array(end)]);
    end

    ylim([0, 75]);
end

