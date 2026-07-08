function [] = plot_energies_with_bands(n, l, m, time_array, ...
        KE_radial_traj_M, KE_theta_traj_M, KE_phi_traj_M, V_traj_M, E_traj_M, use_robust_polar)
    time_ps = time_array / 1e-12;
    % PLOT_ENERGIES_WITH_BANDS  Drop-in replacement for plot_energies.m for
    % runs with M > 1 independent trajectories.
    %
    % Same visual style and analytical-reference yline conventions as
    % plot_energies.m, but each running-average curve is the cross-trajectory
    % mean with a shaded +/- SEM band (SEM = std/sqrt(M)).
    %
    % Optional exception: for the (2,1,0) polar kinetic energy the cross-trajectory
    % mean has formally infinite variance under the analytical Lomax-shape
    % distribution of E_theta (cf. Section 11 of the manuscript), so the
    % curve is the across-trajectory median and the band is the IQR.
    if nargin < 11
        use_robust_polar = true;
    end

    M = size(KE_radial_traj_M, 1);
    hold on;

    % Radial KE -- always Gaussian-tailed -> mean +/- SEM.
    theor_E_r = expected_E_r(n, l);
    if ~isnan(theor_E_r)
        yline(theor_E_r, 'LineStyle', '-', 'Color', [0, 0, 1], 'LineWidth', 0.7, 'HandleVisibility', 'off');
    end
    plot_mean_band(time_ps, KE_radial_traj_M, [0, 0, 1], '-', ...
        '$\left< E_r \right>$');

    % Polar KE -- robust statistics for the heavy-tailed (2,1,0) case.
    yline(expected_E_theta(l, m), 'LineStyle', '-', 'Color', [0, 1, 0], ...
        'LineWidth', 0.7, 'HandleVisibility', 'off');
    if use_robust_polar && n == 2 && l == 1 && m == 0
        % For M = 1 (the manuscript figures) the median/IQR band is degenerate, so
        % drop the "(median, IQR)" qualifier from the legend; keep it on the
        % multi-trajectory dashboards where the band semantics matter.
        if M > 1
            theta_label = '$\left< E_\theta \right>$ (median, IQR)';
        else
            theta_label = '$\left< E_\theta \right>$';
        end
        plot_median_iqr_band(time_ps, KE_theta_traj_M, [0, 1, 0], '--', ...
            theta_label);
    else
        plot_mean_band(time_ps, KE_theta_traj_M, [0, 1, 0], '-', ...
            '$\left< E_\theta \right>$');
    end

    if m ~= 0 % no need to plot trivial energy E_phi=0
        yline(expected_E_phi(l, m), 'LineStyle', '-', 'Color', [0, 1, 0.7], ...
            'LineWidth', 0.7, 'HandleVisibility', 'off');
        plot_mean_band(time_ps, KE_phi_traj_M, [0, 1, 0.7], '-', ...
            '$\left< E_\phi \right>$');
    end

    yline(expected_V(n), 'LineStyle', '-', 'Color', [1, 0.5, 0], ...
        'LineWidth', 0.5, 'HandleVisibility', 'off');
    plot_mean_band(time_ps, V_traj_M, [1, 0.5, 0], '-', ...
        '$\left< V \right>$');

    yline(expected_E(n), 'LineStyle', '-', 'Color', [1, 0, 0], ...
        'LineWidth', 0.7, 'HandleVisibility', 'off');
    plot_mean_band(time_ps, E_traj_M, [1, 0, 0], '-', ...
        '$\left< E \right>$');

    hold off;
    ylabel('Average Energies (J)', 'FontSize', 14);
    xlabel('Time (ps)', 'FontSize', 14);
    legend('Interpreter', 'latex', 'FontSize', 11);
    set(gca, 'FontSize', 12);

    if length(time_ps) > 1
        xlim([0, time_ps(end)]);
    end

    % Match plot_energies.m for the manuscript's 2p_0 and 2s_0 figures.
    if (n == 2 && l == 1 && m == 0) || (n == 2 && l == 0 && m == 0)
        yticks(-15e-19:3e-19:15e-19);
        ylim([-15e-19, 15e-19]);
    end

    % Annotate the ensemble size so the reader knows the band semantics, but only
    % for multi-trajectory runs (M > 1, the cutoff-scan dashboards). For the
    % single-trajectory manuscript figures (M = 1) there are no bands, so the
    % annotation is omitted -- state M = 1 in the figure caption instead.
    if M > 1
        text(0.02, 0.02, sprintf('M = %d trajectories', M), ...
            'Units', 'normalized', 'FontSize', 10, ...
            'BackgroundColor', [1 1 1 0.7], 'Margin', 2);
    end
end

function plot_mean_band(t, X_M, color, line_style, display_name)
    t = t(:).';
    mu  = mean(X_M, 1);
    sem = std(X_M, 0, 1) / sqrt(size(X_M, 1));
    band_x = [t, fliplr(t)];
    band_y = [mu + sem, fliplr(mu - sem)];
    patch(band_x, band_y, color, ...
        'FaceAlpha', 0.18, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    plot(t, mu, 'LineStyle', line_style, 'Color', color, ...
        'LineWidth', 2, 'DisplayName', display_name);
end

function plot_median_iqr_band(t, X_M, color, line_style, display_name)
    t = t(:).';
    med = median(X_M, 1);
    q1  = quantile(X_M, 0.25, 1);
    q3  = quantile(X_M, 0.75, 1);
    band_x = [t, fliplr(t)];
    band_y = [q3, fliplr(q1)];
    patch(band_x, band_y, color, ...
        'FaceAlpha', 0.18, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    plot(t, med, 'LineStyle', line_style, 'Color', color, ...
        'LineWidth', 2, 'DisplayName', display_name);
end

% --- Analytical references (kept in sync with plot_energies.m) -------------
function [exp_E_r] = expected_E_r(n, l)
    [hbar, m_e, a_0, ~, ~, ~] = constants();
    mean_inv_r2 = 1 / (12 * a_0^2);
    exp_E_r = hbar^2 / (2 * m_e) * (1 / (a_0^2 * n^2) - mean_inv_r2 * l * (l + 1));
end

function [exp_E_theta] = expected_E_theta(l, m)
    [hbar, m_e, a_0, ~, ~, ~] = constants();
    mean_inv_r2 = 1 / (12 * a_0^2);
    m_abs = abs(m);
    exp_E_theta = hbar^2 / (2 * m_e) * mean_inv_r2 * (2 * l^2 - 2 * l * (m_abs - 1) - m_abs) / 2;
end

function [exp_E_phi] = expected_E_phi(l, m)
    [hbar, m_e, a_0, ~, ~, ~] = constants();
    mean_inv_r2 = 1 / (12 * a_0^2);
    exp_E_phi = hbar^2 / (2 * m_e) * mean_inv_r2 * (2 * l + 1) * abs(m) / 2;
end

function [exp_V] = expected_V(n)
    exp_V = 2 * expected_E(n);
end

function [exp_E] = expected_E(n)
    [hbar, m_e, a_0, ~, ~, ~] = constants();
    exp_E = - hbar^2 / (2 * m_e) / (a_0^2 * n^2);
end
