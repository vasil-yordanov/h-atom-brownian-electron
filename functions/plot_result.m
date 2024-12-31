function plot_result(v, S, simTime, n_steps, dt, frame_index, ...
    hist_counts_r, hist_bins_r, ...
    hist_counts_theta, hist_bins_theta, ...
    hist_counts_phi, hist_bins_phi, ...
    hist_bins_dS, hist_counts_dS, ...
    r_values, P_analytic_r, ...
    theta_values, P_analytic_theta, ...
    phi_values, P_analytic_phi ...
)
        a0 = 5.29177e-11; 
        fprintf('%. ->  %.5e\r', simTime, S/simTime);

        % Normalize the Histogram Counts
        hist_counts_r_norm = normalize_histogram(hist_counts_r, hist_bins_r);
        hist_counts_theta_norm = normalize_histogram(hist_counts_theta, hist_bins_theta);
        hist_counts_phi_norm = normalize_histogram(hist_counts_phi, hist_bins_phi);

        % Store the Action and Time at the Current Frame
        S_traj(frame_index) = S;

        time_array(frame_index) = simTime;

        % Compute Deviations for R, theta, and phi Distributions
        % Calculate bin centers for R histogram
        r_bin_centers = hist_bins_r(1:end-1) + diff(hist_bins_r)/2;
        P_analytic_r_bin = interp1(r_values, P_analytic_r, r_bin_centers, 'linear', 'extrap');

        deviation_R = sqrt(sum((hist_counts_r_norm - P_analytic_r_bin) .^2) * mean(diff(r_bin_centers)));

        % Calculate bin centers for theta histogram
        theta_bin_centers = hist_bins_theta(1:end-1) + diff(hist_bins_theta)/2;
        P_analytic_theta_bin = interp1(theta_values, P_analytic_theta, theta_bin_centers, 'linear', 'extrap');

        deviation_theta = sqrt(sum((hist_counts_theta_norm - P_analytic_theta_bin) .^2) * mean(diff(theta_bin_centers)));

        % Calculate bin centers for phi histogram
        phi_bin_centers = hist_bins_phi(1:end-1) + diff(hist_bins_phi)/2;
        P_analytic_phi_bin = interp1(phi_values, P_analytic_phi, phi_bin_centers, 'linear', 'extrap');

        deviation_phi = sqrt(sum((hist_counts_phi_norm - P_analytic_phi_bin) .^2) * mean(diff(phi_bin_centers)));

        % Store Deviations
        deviation_array_R(frame_index) = deviation_R;
        deviation_array_theta(frame_index) = deviation_theta;
        deviation_array_phi(frame_index) = deviation_phi;

        % Create a Figure with Six Subplots
        figure(100);  % Use a specific figure number
        clf;  % Clear the figure

        % Set Figure Size to be Optimized for Screen
        set(gcf, 'Position', [50, 50, 1600, 900]);  % [left, bottom, width, height]
        % After creating subplots, adjust the layout

        %% --- Top Row: R, Theta, and Phi Histograms ---
        % --- Left Subplot: Radial Distribution ---
        subplot(2,4,1);
        bar(r_bin_centers / a0, hist_counts_r_norm, 'FaceAlpha', 0.6, 'EdgeColor', 'none');
        hold on;
        plot(r_values / a0, P_analytic_r, 'r', 'LineWidth', 2);
        hold off;
        xlabel('r / a_0');
        ylabel('Normalized Radial Probability Density');
        title('Radial (R) Distribution');
        legend('Simulation', 'Analytic Solution');
        set(gca, 'FontSize', 12);

        % --- Middle Subplot: Theta Distribution ---
        subplot(2,4,2);
        bar(theta_bin_centers, hist_counts_theta_norm, 'FaceAlpha', 0.6, 'EdgeColor', 'none');
        hold on;
        plot(theta_values, P_analytic_theta, 'r', 'LineWidth', 2);
        hold off;
        xlabel('\theta (radians)');
        ylabel('Normalized Polar Probability Density');
        title('Polar (\theta) Distribution');
        legend('Simulation', 'Analytic Solution');
        set(gca, 'FontSize', 12);

        % --- Right Subplot: Azimuthal Distribution ---
        subplot(2,4,3);
        bar(phi_bin_centers, hist_counts_phi_norm, 'FaceAlpha', 0.6, 'EdgeColor', 'none');
        hold on;
        plot(phi_values, P_analytic_phi, 'r', 'LineWidth', 2);
        hold off;
        xlabel('\phi (radians)');
        ylabel('Normalized Azimuthal Probability Density');
        title('Azimuthal (\phi) Distribution');
        legend('Simulation', 'Analytic Solution');
        set(gca, 'FontSize', 12);

        %% --- Bottom Row: Deviations and Action over Time ---
        % --- Combined Deviation Plot with Dual Y-Axes ---
        subplot(2,4,4);
        yyaxis left
        plot(time_array(1:frame_index), deviation_array_R(1:frame_index), 'b', 'DisplayName', 'Radial (R) Deviation');
        ylabel('Deviation R');
        hold on;

        yyaxis right
        plot(time_array(1:frame_index), deviation_array_theta(1:frame_index), 'g', 'DisplayName', 'Polar (\theta) Deviation');
        plot(time_array(1:frame_index), deviation_array_phi(1:frame_index), 'm', 'DisplayName', 'Azimuthal (\phi) Deviation');
        ylabel('Deviation \theta and \phi');

        hold off;
        set(gca, 'FontSize', 12);
        xlabel('Time (s)');
        title('Deviation of Distributions over Time');
        legend('Location', 'best');
        xlim([0, n_steps * dt]);

        % Adjust Y-axis Limits
        yyaxis left
        ylim_left = [0, max(deviation_array_R(1:frame_index)) * 1.1];
        ylim(ylim_left);

        yyaxis right
        ylim_right = [0, max([deviation_array_theta(1:frame_index), deviation_array_phi(1:frame_index)]) * 1.1];
        ylim(ylim_right);

        % --- Right Subplot: Action per Time over Time ---
        subplot(2,4,5);
        % Plot the action divided by time as a trace over time
        action_over_time = S_traj(1:frame_index) ./ time_array(1:frame_index);
        % Handle division by zero at time zero
        action_over_time(time_array(1:frame_index) == 0) = 0;  % Define action/time at t=0 as zero
        plot(time_array(1:frame_index), action_over_time, 'b');
        set(gca, 'FontSize', 12);
        xlabel('Time (s)');
        ylabel('Action per Time (J)');
        title('J/t (t)');
        xlim([0, n_steps * dt]);
        % Adjust y-axis limits dynamically

        % --- Sixth Subplot: Log-Log Plot of dS ---
        subplot(2,4,6);

        % Remove zero or negative values from hist_counts_dS
        positive_dS_bins = hist_bins_dS(1:end-1);  % Use the correct bin centers (exclude the last edge)
        positive_dS_counts = hist_counts_dS(hist_counts_dS > 0);  % Remove empty bins (zero counts)

        % Ensure that we are only working with bins that have counts
        positive_dS_bins = positive_dS_bins(hist_counts_dS > 0);  % Match the size to the filtered counts

        % Filter out zero or negative values (to avoid NaNs in log)
        valid_idx = (positive_dS_bins > 0) & (positive_dS_counts > 0);
        positive_dS_bins = positive_dS_bins(valid_idx);
        positive_dS_counts = positive_dS_counts(valid_idx);

        % Debugging: print out min values to ensure no invalid data
        fprintf('Min positive_dS_bins: %.2e\n', min(positive_dS_bins));
        fprintf('Min positive_dS_counts: %.2e\n', min(positive_dS_counts));

        if ~isempty(positive_dS_counts) && length(positive_dS_bins) == length(positive_dS_counts)
            % Log-log plot of dS
            loglog(positive_dS_bins, positive_dS_counts, 'o');
            hold on;

            % Fit a power law y = a*x^alpha if there are enough data points
            if length(positive_dS_bins) > 1
                % Polyfit in log-log space
                coeffs = polyfit(log10(positive_dS_bins), log10(positive_dS_counts), 1);
                alpha = -coeffs(1);
                a_fit = 10^coeffs(2);
                alpha_values(frame_index) = alpha;

                % Calculate the fitted line over the same range
                fitted_line = a_fit * positive_dS_bins.^(-alpha);

                % Plot the fitted line
                loglog(positive_dS_bins, fitted_line, 'r--', 'LineWidth', 2);

                % Output fit parameters for debugging
                fprintf('Fit Parameters: a = %.2f, alpha = %.2f\n', a_fit, alpha);

                % Check if fitted_line values make sense
                fprintf('Fitted line min: %.2e, max: %.2e\n', min(fitted_line), max(fitted_line));
                 % Add labels and title
                xlabel('dS');
                ylabel('Frequency');
                title('Log-Log Plot of dS');
                legend('Data', sprintf('Fit: y = %.2f x^{%.2f}', a_fit, -alpha), 'Location', 'best');

                % Force the plot to use an appropriate y-axis range
                ylim([min(positive_dS_counts)/10, max(positive_dS_counts)*10]);
                xlim([min(positive_dS_bins)/10, max(positive_dS_bins)*10]);

            end


            hold off;
        end

        subplot(2,4,7);  % Adjust the subplot index to be the 7th plot
        plot(time_array(1:frame_index), alpha_values(1:frame_index), 'LineWidth', 2);  % Time in ps
        xlabel('Time (s)');
        ylabel('\alpha');
        title('Evolution of \alpha Over Time');
        ylim([1, 4]);  % Set y-axis limits based on alpha range
        xlim([0, n_steps * dt]);

        % --- Add Simulation Time Text ---
        annotation_text = sprintf('t = %.2e s', simTime);
        annotation('textbox', [0.02, 0.95, 0.1, 0.05], 'String', annotation_text, ...
            'Units', 'normalized', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', ...
            'FontSize', 12, 'FontWeight', 'bold', 'EdgeColor', 'none');

        % Capture the Frame
        frame = getframe(gcf);
        writeVideo(v, frame);
end

