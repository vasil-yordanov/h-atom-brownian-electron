function [] = analyse_cutoff_scan(state_label, export_dir, data_root, mode, mix_time)
    % ANALYSE_CUTOFF_SCAN  Kinetic-energy cutoff scan post-processing.
    %
    % Usage:
    %   analyse_cutoff_scan('2s0');                 % default: manuscript Fig. 11 only
    %   analyse_cutoff_scan('2s0', 'figures');
    %   analyse_cutoff_scan('2s0', 'figures', 'data', 'extended');
    %   analyse_cutoff_scan('2s0', 'figures', 'data', 'default', true);  % mix total times
    %
    % Loads every scan .mat file under data_root (all seeds) and checks whether
    % the relevant kinetic-energy first moment is stable under changes in v_max.
    %
    % mode:
    %   'default'  (default) -- print the convergence tables and export only
    %              the single-panel manuscript figure (Fig. 11).
    %   'extended' -- additionally export a 6-panel diagnostic figure:
    %              (a) mean +/- bootstrap SE, (b) median +/- bootstrap SE,
    %              (c) all per-trajectory points, (d) nodal crossing rate,
    %              (e) running-mean error, (f) running trajectory scatter.
    %              The mean/median panels deliberately show both mean/SEM and
    %              median/IQR views because the nodal kinetic energy is
    %              heavy-tailed.
    %
    % mix_time (default false) -- DIAGNOSTIC. When false the scan is expected to
    %   share a single total simulated time and curves are grouped/labelled by
    %   the integration step Delta t only. When true, runs with different total
    %   times T are allowed on the same figure: curves are grouped by (Delta t, T)
    %   so each connected line is still a single total time, the legend shows
    %   both Delta t and T, and the output files get a '_mixtime' suffix so the
    %   single-time manuscript figure is not overwritten. The per-trajectory mean
    %   <T> is a per-step time average, so it is comparable across total times;
    %   mixing is therefore meaningful for the mean/median/scatter panels.

    if nargin < 2
        export_dir = 'figures';
    end
    if nargin < 3
        data_root = 'data';
    end
    if nargin < 4
        mode = 'default';
    end
    if nargin < 5 || isempty(mix_time)
        mix_time = false;
    end
    if ischar(mix_time) || isstring(mix_time)
        mix_time = any(strcmpi(char(mix_time), ...
            {'mix', 'mixtime', 'mix_time', '--mix-max-time', 'true'}));
    end
    if ~(strcmp(mode, 'default') || strcmp(mode, 'extended'))
        error('mode must be either ''default'' or ''extended''.');
    end
    extended = strcmp(mode, 'extended');
    name_suffix = '';
    if mix_time
        name_suffix = '_mixtime';
    end
    if ~exist(export_dir, 'dir')
        mkdir(export_dir);
    end

    if ~(strcmp(state_label, '2p0') || strcmp(state_label, '2s0'))
        error('state_label must be either ''2p0'' or ''2s0''.');
    end

    % Load every scan run for this state (all seeds). Folder/file convention,
    % with the seed as the last token:
    %   2s0_scan_vmax_3p0c_seed_17000
    %   2s0_scan_vmax_3p0c_dt_1zs_seed_31000
    files = dir(fullfile(data_root, [state_label '_scan_vmax_*'], ...
                         [state_label '_scan_vmax_*.mat']));
    if isempty(files)
        error('No scan .mat files found for state %s in %s.', state_label, data_root);
    end

    % One curve per dt group, one point per (dt, v_max) cell: keep the run
    % with the LONGEST total simulated time in each cell (best statistics)
    % and ignore any further seed replicates.
    % Untagged runs use the default dt = 10 zs (1e-20 s).
    pre_v = zeros(numel(files), 1);
    pre_dt = zeros(numel(files), 1);
    pre_T = zeros(numel(files), 1);
    for k = 1:numel(files)
        [pre_v(k), ~, pre_dt(k)] = parse_vmax_from_filename(state_label, files(k).name);
        try
            q = load(fullfile(files(k).folder, files(k).name), 'n_steps', 'dt');
            pre_T(k) = q.n_steps * q.dt;
        catch
            pre_T(k) = 0;
        end
    end
    [~, ord] = sort(pre_T, 'descend');
    files = files(ord);
    pre_v = pre_v(ord);
    pre_dt = pre_dt(ord);
    pre_T = pre_T(ord);
    if mix_time
        % Keep best seed per (dt, v_max, T-bucket); different total times
        % are distinct curve points and must not be merged.
        t_bucket = round(pre_T / 1e-13);
        [~, keep_first] = unique([pre_dt, pre_v, t_bucket], 'rows', 'stable');
    else
        [~, keep_first] = unique([pre_dt, pre_v], 'rows', 'stable');
    end
    files = files(keep_first);

    n_files = length(files);
    v_over_c = zeros(n_files, 1);
    T_mean = zeros(n_files, 1);
    T_sem = zeros(n_files, 1);
    T_median = zeros(n_files, 1);
    T_mean_se = zeros(n_files, 1);    % bootstrap SE of the sample mean
    T_median_se = zeros(n_files, 1);  % bootstrap SE of the sample median
    T_q1 = zeros(n_files, 1);
    T_q3 = zeros(n_files, 1);
    T_rel_mean = zeros(n_files, 1);
    T_rel_median = zeros(n_files, 1);
    cutoff_engagement_mean = zeros(n_files, 1);
    N_cross_mean = NaN(n_files, 1);
    N_cross_std = NaN(n_files, 1);
    N_cross_total = NaN(n_files, 1);
    seed_id = NaN(n_files, 1);
    dt_zs = NaN(n_files, 1);   % integration step in zeptoseconds (default 10 = 1e-20 s)
    T_max_per_traj = zeros(n_files, 1);
    T_min_per_traj = zeros(n_files, 1);
    M_count = zeros(n_files, 1);
    T_per_traj_cell = cell(n_files, 1);   % full per-trajectory vectors (scatter panel)
    total_time = NaN(n_files, 1);         % total simulated time T = time_arr(end)

    % Bootstrap settings for the seed-to-seed estimator-uncertainty
    % error bars on panels (a) and (b). The heavy Lomax tail of T_r
    % (2s0) and T_theta (2p0) means the classical SEM/IQR are not
    % faithful uncertainty quantifiers. See README_analyse_cutoff_scan.md.
    n_boot = 2000;
    boot_rng_seed = 17;

    time_arr_ref = [];
    T_mean_time = [];
    T_sem_time = [];

    T_analytical = kinetic_energy_reference_value(state_label);

    for k = 1:n_files
        mat_name = files(k).name;
        mat_path = fullfile(files(k).folder, mat_name);
        [v_over_c(k), seed_id(k), dt_zs(k)] = parse_vmax_from_filename(state_label, mat_name);

        S = load(mat_path, 'S_KE_radial_M', 'S_KE_theta_M', ...
                           'KE_radial_traj_M', 'KE_theta_traj_M', ...
                           'cutoff_engagement_r_per_traj', ...
                           'cutoff_engagement_theta_per_traj', ...
                           'N_cross_per_traj', 'N_cross_total', ...
                           'n_steps', 'time_arr');
        require_field(S, 'n_steps', mat_name);

        if strcmp(state_label, '2p0')
            require_field(S, 'S_KE_theta_M', mat_name);
            require_field(S, 'cutoff_engagement_theta_per_traj', mat_name);
            T_per_traj = S.S_KE_theta_M / S.n_steps;
            engagement = S.cutoff_engagement_theta_per_traj;
            T_time = optional_time_series(S, 'KE_theta_traj_M');
        else
            require_field(S, 'S_KE_radial_M', mat_name);
            require_field(S, 'cutoff_engagement_r_per_traj', mat_name);
            T_per_traj = S.S_KE_radial_M / S.n_steps;
            engagement = S.cutoff_engagement_r_per_traj;
            T_time = optional_time_series(S, 'KE_radial_traj_M');
        end

        T_mean(k) = mean(T_per_traj);
        T_sem(k) = std(T_per_traj) / sqrt(length(T_per_traj));
        T_median(k) = median(T_per_traj);
        T_q1(k) = quantile(T_per_traj, 0.25);
        T_q3(k) = quantile(T_per_traj, 0.75);
        % Honest seed-to-seed estimator uncertainty for both panels:
        T_mean_se(k)   = bootstrap_mean_se(T_per_traj, n_boot, boot_rng_seed + k);
        T_median_se(k) = bootstrap_median_se(T_per_traj, n_boot, boot_rng_seed + k);
        T_rel_mean(k) = abs(T_mean(k) - T_analytical) / T_analytical;
        T_rel_median(k) = abs(T_median(k) - T_analytical) / T_analytical;
        cutoff_engagement_mean(k) = mean(engagement);
        T_max_per_traj(k) = max(T_per_traj);
        T_min_per_traj(k) = min(T_per_traj);
        M_count(k) = length(T_per_traj);
        T_per_traj_cell{k} = T_per_traj(:);
        if isfield(S, 'time_arr') && ~isempty(S.time_arr)
            total_time(k) = S.time_arr(end);
        end
        if isfield(S, 'N_cross_per_traj')
            N_cross_mean(k) = mean(S.N_cross_per_traj);
            N_cross_std(k) = std(S.N_cross_per_traj);
            N_cross_total(k) = sum(S.N_cross_per_traj);
        elseif isfield(S, 'N_cross_total')
            N_cross_total(k) = S.N_cross_total;
        end

        if ~isempty(T_time)
            require_field(S, 'time_arr', mat_name);
            if isempty(time_arr_ref)
                time_arr_ref = S.time_arr;
                T_mean_time = NaN(n_files, length(time_arr_ref));
                T_sem_time = NaN(n_files, length(time_arr_ref));
            else
                % Compare time grids up to floating-point round-off. dt sweeps
                % hold the total time fixed but set n_steps = round(total/dt),
                % so grids that are physically identical can differ by ~1e-14
                % relative. Runs with a genuinely different total simulated
                % time are kept for the (time-independent) per-trajectory
                % averages, but their running time series is skipped on the
                % time-series diagnostic panels.
                grid_scale = max(abs(time_arr_ref(:)));
                if length(S.time_arr) ~= length(time_arr_ref) || ...
                        any(abs(S.time_arr(:) - time_arr_ref(:)) > 1e-9 * grid_scale)
                    warning('Time grid differs in %s; skipping its running time series.', ...
                            mat_name);
                    T_time = [];
                end
            end
        end
        if ~isempty(T_time)
            T_mean_time(k, :) = mean(T_time, 1); %#ok<AGROW>
            T_sem_time(k, :) = std(T_time, 0, 1) / sqrt(size(T_time, 1)); %#ok<AGROW>
        end
    end

    [v_over_c, idx] = sort(v_over_c);
    T_mean = T_mean(idx);
    T_sem = T_sem(idx);
    T_median = T_median(idx);
    T_mean_se = T_mean_se(idx);
    T_median_se = T_median_se(idx);
    T_q1 = T_q1(idx);
    T_q3 = T_q3(idx);
    T_rel_mean = T_rel_mean(idx);
    T_rel_median = T_rel_median(idx);
    cutoff_engagement_mean = cutoff_engagement_mean(idx);
    N_cross_mean = N_cross_mean(idx);
    N_cross_std = N_cross_std(idx);
    N_cross_total = N_cross_total(idx);
    seed_id = seed_id(idx);
    dt_zs = dt_zs(idx);
    T_max_per_traj = T_max_per_traj(idx);
    T_min_per_traj = T_min_per_traj(idx);
    M_count = M_count(idx);
    T_per_traj_cell = T_per_traj_cell(idx);
    total_time = total_time(idx);
    if ~isempty(T_mean_time)
        T_mean_time = T_mean_time(idx, :);
        T_sem_time = T_sem_time(idx, :);
    end

    % Curve grouping for the figures: by integration step Delta t alone
    % (default), or by (Delta t, total time T) when mix_time is on so that
    % runs of different total simulated time stay on separate, correctly
    % labelled curves. group_id(k) is the group of point k; group_labels{g}
    % is its legend string.
    [group_id, group_labels] = make_curve_groups(dt_zs, total_time, mix_time);

    % Derived diagnostic quantities
    alpha_fs    = 7.2973525693e-3;           % fine-structure constant
    d_cut_a0    = alpha_fs ./ v_over_c;      % nodal-layer width in units of a_0
    bias_mean   = (T_mean   - T_analytical) / T_analytical;  % signed
    bias_median = (T_median - T_analytical) / T_analytical;  % signed
    z_mean      = bias_mean   ./ max(T_mean_se   / T_analytical, eps);
    z_median    = bias_median ./ max(T_median_se / T_analytical, eps);

    fprintf('\n=== Kinetic cutoff scan: %s  |  %s ===\n', ...
            state_label, data_root);
    fprintf('Analytical T = %.6e J\n', T_analytical);
    fprintf('bias = (T - T_an)/T_an  (signed).  z = bias / bootSE  (|z|<2 consistent with analytical).\n');
    fprintf('d_cut/a0 = alpha_fs / (v_max/c)  (nodal-layer width in Bohr radii).\n\n');

    % Table 1: per-run diagnostics
    hdr = '  %7s %6s %7s %4s %10s %10s %8s %10s %8s %8s %10s %10s %10s %10s';
    row = '  %7s %6.3g %7g %4d %10.3e %10.3e %8.2f %10.3e %8.2f %8.3e %10.3e %10.3e %10.2e %10.2e';
    fprintf(hdr, 'seed', 'vmax/c', 'dt_zs', 'M', ...
            'bias_mean', 'SEm/T_an', 'z_mean', ...
            'bias_med',  'SEmed/T_an', 'z_med', ...
            'd_cut/a0', 'engagement', 'Tmax/T_an', 'Tmin/T_an');
    fprintf('\n');
    for k = 1:length(v_over_c)
        if isnan(seed_id(k))
            slabel = 'common';
        else
            slabel = sprintf('s%d', seed_id(k));
        end
        fprintf(row, slabel, v_over_c(k), dt_zs(k), M_count(k), ...
                bias_mean(k),   T_mean_se(k)   / T_analytical, z_mean(k), ...
                bias_median(k), T_median_se(k) / T_analytical, z_median(k), ...
                d_cut_a0(k), cutoff_engagement_mean(k), ...
                T_max_per_traj(k) / T_analytical, T_min_per_traj(k) / T_analytical);
        fprintf('\n');
    end
    fprintf('\n');

    % Table 2: cross-seed spread grouped by v_max/c
    unique_v = unique(v_over_c);
    fprintf('=== Cross-seed spread by v_max/c (median values, in units of T_an) ===\n');
    fprintf('  spread/SE > 2 suggests the cutoff is too loose for M trajectories to self-average.\n');
    fprintf('  %8s %4s %10s %10s %10s %10s %10s\n', ...
            'vmax/c', 'N', 'med_min', 'med_max', 'spread', 'mean_SEmed', 'spread/SE');
    for u = 1:length(unique_v)
        rows = find(v_over_c == unique_v(u));
        meds = T_median(rows) / T_analytical;
        ses  = T_median_se(rows) / T_analytical;
        sprd = max(meds) - min(meds);
        mse  = mean(ses);
        fprintf('  %8.3g %4d %10.4f %10.4f %10.4f %10.4f %10.3g\n', ...
                unique_v(u), length(rows), min(meds), max(meds), sprd, mse, sprd/max(mse,eps));
    end
    fprintf('=======================================\n');

    if ~extended
        % -------------------------------------------------------------
        % Default deliverable: the manuscript figure (Physica Scripta
        % style, Fig. 11) -- a single clean panel showing the mean only,
        % with one connected error-bar curve per dt group.
        % -------------------------------------------------------------
        export_manuscript_mean_figure(state_label, export_dir, name_suffix, ...
            v_over_c, T_mean, T_mean_se, group_id, group_labels, T_analytical);
        return;
    end

    % -----------------------------------------------------------------
    % Extended diagnostic figure: 6 panels in a 3x2 grid.
    %   (a) mean   +/- bootstrap SE   (whiskers: classical SEM)
    %   (b) median +/- bootstrap SE   (whiskers: IQR)
    %   (c) all per-trajectory points (heavy-tail / seed spread)
    %   (d) nodal crossing rate per trajectory
    %   (e) running-mean error vs time (or final error vs v_max)
    %   (f) running trajectory scatter vs time (or final SEM vs v_max)
    % -----------------------------------------------------------------
    fig = figure('Position', [50, 50, 1500, 1300]);
    set(fig, 'DefaultAxesFontSize', 13, ...
             'DefaultTextFontSize', 13, ...
             'DefaultLegendFontSize', 12);

    colors = line_colors(length(v_over_c));

    subplot(3, 2, 1);
    % Primary error bars: bootstrap SE of the sample mean. The classical
    % SEM uses the empirical sample variance, which under the Lomax tail
    % (T_r for 2s0; T_theta for 2p0) is bounded only by the cutoff and
    % therefore understates true seed-to-seed uncertainty as v_max/c grows.
    plot_dt_grouped_errorbar(v_over_c, T_mean, T_mean_se, T_mean_se, group_id, group_labels, 'o');
    % Light whiskers: the classical SEM, kept as a descriptive overlay.
    plot_light_whiskers(v_over_c, T_mean, T_sem, T_sem);
    set(gca, 'XScale', 'log');
    xlabel('$v_{\max}/c$', 'Interpreter', 'latex');
    ylabel('$\left< T \right>$ (J)', 'Interpreter', 'latex');
    title('(a) Mean $\pm$ bootstrap SE (whiskers: classical SEM)', ...
          'Interpreter', 'latex');
    grid on;
    hold on;
    draw_horizontal_line(T_analytical, '-.', [0.60, 0.00, 0.00]);
    legend('show', 'Interpreter', 'latex', 'Location', 'best', 'Box', 'off');
    hold off;

    subplot(3, 2, 2);
    % Primary error bars: bootstrap SE of the sample median. IQR is a
    % within-sample spread statistic (converges to a finite population
    % IQR as M -> Inf) and does NOT quantify the uncertainty of the
    % median estimator across different draws of M trajectories.
    plot_dt_grouped_errorbar(v_over_c, T_median, T_median_se, T_median_se, group_id, group_labels, 's');
    % Light whiskers: the IQR, kept as a descriptive overlay of the
    % within-sample spread.
    plot_light_whiskers(v_over_c, T_median, T_median - T_q1, T_q3 - T_median);
    set(gca, 'XScale', 'log');
    xlabel('$v_{\max}/c$', 'Interpreter', 'latex');
    ylabel('$\mathrm{median}\left(\left< T \right>_j\right)$ (J)', 'Interpreter', 'latex');
    title('(b) Median $\pm$ bootstrap SE (whiskers: IQR)', ...
          'Interpreter', 'latex');
    grid on;
    hold on;
    draw_horizontal_line(T_analytical, '-.', [0.60, 0.00, 0.00]);
    legend('show', 'Interpreter', 'latex', 'Location', 'best', 'Box', 'off');
    hold off;

    subplot(3, 2, 3);
    % Every trajectory's time-averaged kinetic energy as one point, grouped
    % by dt (colour) at x = v_max/c, with per-file mean (diamond) and median
    % (tick) overlaid. Exposes the heavy (Lomax) tail and seed-to-seed
    % spread directly, instead of summarising them with a single bar.
    plot_trajectory_points_panel(v_over_c, T_per_traj_cell, group_id, group_labels, T_analytical);
    title('(c) All per-trajectory $\left< T \right>$ (diamond: mean, tick: median)', ...
          'Interpreter', 'latex');

    subplot(3, 2, 4);
    if any(~isnan(N_cross_mean))
        plot_dt_grouped_errorbar(v_over_c, N_cross_mean, N_cross_std, N_cross_std, group_id, group_labels, 'o');
        set(gca, 'XScale', 'log');
        if all(N_cross_mean(~isnan(N_cross_mean)) > 0)
            set(gca, 'YScale', 'log');
        end
        xlabel('$v_{\max}/c$', 'Interpreter', 'latex');
        ylabel('$N_{\rm cross}$ per trajectory', 'Interpreter', 'latex');
        title('(d) Nodal crossing rate', 'Interpreter', 'latex');
        grid on;
        hold on;
        draw_horizontal_line(100, '--', [0.60, 0.00, 0.00]);
        legend('show', 'Interpreter', 'latex', 'Location', 'best', 'Box', 'off');
        hold off;
    else
        text(0.5, 0.5, 'no crossing data', 'Units', 'normalized', ...
             'HorizontalAlignment', 'center', 'Interpreter', 'latex');
        title('(d) Nodal crossing rate', 'Interpreter', 'latex');
        axis off;
    end

    subplot(3, 2, 5);
    hold on;
    if isempty(T_mean_time)
        plot_dt_grouped_line(v_over_c, T_rel_mean, group_id, group_labels, 'o');
        set(gca, 'XScale', 'log');
        xlabel('$v_{\max}/c$', 'Interpreter', 'latex');
        ylabel('final relative error', 'Interpreter', 'latex');
    else
        t_ps = time_arr_ref / 1e-12;
        for k = 1:length(v_over_c)
            line_style = line_style_for_seed(seed_id(k));
            plot(t_ps, abs(T_mean_time(k, :) - T_analytical) / T_analytical, ...
                 line_style, 'Color', colors(k, :), 'LineWidth', 1.2);
        end
        xlabel('$t$ (ps)', 'Interpreter', 'latex');
        ylabel('$|\overline{T}(t)-T_{\mathrm{analytic}}|/T_{\mathrm{analytic}}$', 'Interpreter', 'latex');
    end
    title('(e) Running mean error', 'Interpreter', 'latex');
    grid on;
    hold off;

    subplot(3, 2, 6);
    hold on;
    if isempty(T_sem_time)
        plot_dt_grouped_line(v_over_c, T_sem / T_analytical, group_id, group_labels, '^');
        set(gca, 'XScale', 'log');
        xlabel('$v_{\max}/c$', 'Interpreter', 'latex');
        legend('show', 'Interpreter', 'latex', ...
               'Location', 'northeastoutside', 'Box', 'off');
    else
        t_ps = time_arr_ref / 1e-12;
        for k = 1:length(v_over_c)
            line_style = line_style_for_seed(seed_id(k));
            plot(t_ps, T_sem_time(k, :) / T_analytical, ...
                 line_style, 'Color', colors(k, :), 'LineWidth', 1.2, ...
                 'HandleVisibility', 'off');
        end
        xlabel('$t$ (ps)', 'Interpreter', 'latex');
    end
    ylabel('$\mathrm{SEM}\left(\left< T \right>_j\right)/T_{\mathrm{analytic}}$', 'Interpreter', 'latex');
    title('(f) Running trajectory scatter', 'Interpreter', 'latex');
    grid on;
    hold off;

    title_str = sprintf('Kinetic-energy cutoff scan (%s)', state_label);
    try
        sgtitle(title_str, 'FontSize', 16, 'FontWeight', 'bold');
    catch
        annotation('textbox', [0, 0.96, 1, 0.04], 'String', title_str, ...
                   'HorizontalAlignment', 'center', 'EdgeColor', 'none', ...
                   'FontSize', 16, 'FontWeight', 'bold');
    end

    pdf_name = [state_label '_cutoff_scan_kinetic' name_suffix '.pdf'];
    exportgraphics(fig, fullfile(export_dir, pdf_name), 'ContentType', 'vector');
    fprintf('%s written\n', fullfile(export_dir, pdf_name));
end

function export_manuscript_mean_figure(state_label, export_dir, name_suffix, ...
                                       v_over_c, T_mean, T_mean_se, group_id, group_labels, T_analytical)
    % Clean single-panel mean(<T>) vs v_max/c figure for the manuscript. Each of
    % the two states is included at width=0.49\textwidth (= 204 pt). exportgraphics
    % crops the vector PDF to its content, so the figure is built deliberately LARGE
    % (page well above 204 pt) and is then scaled DOWN by LaTeX, exactly like the
    % distribution/energy figures. Down-scaling a vector PDF stays crisp; the
    % earlier small builds were blown UP ~2x into the slot, which softened the text.
    % Font sizes are picked so the down-scaled labels (~12 pt) match those figures.
    % LaTeX only on labels/legend (default tick numbers) keeps it crisp; forcing the
    % LaTeX tick interpreter would embed Computer Modern as unhinted font subsets
    % that look soft. One connected curve per dt group.

    fig = figure('Units', 'centimeters', 'Position', [2, 2, 20, 15.5], ...
                 'Color', 'w');
    ax = axes('Parent', fig);
    set(ax, 'FontSize', 12, 'Box', 'on', 'LineWidth', 1.0);
    hold(ax, 'on');

    % Print-safe palette and distinct markers, one per curve group.
    palette = [0.00, 0.20, 0.55;    % blue
               0.10, 0.45, 0.20;    % green
               0.70, 0.30, 0.00;    % orange-brown
               0.45, 0.10, 0.55;    % purple
               0.10, 0.55, 0.60];   % teal (5th group, e.g. 1 zs)
    markers = {'o', 's', '^', 'd', 'v'};

    for g = 1:numel(group_labels)
        gi = (group_id == g);
        if ~any(gi), continue; end
        [vs, order] = sort(v_over_c(gi));
        yy = T_mean(gi);    yy = yy(order);
        ee = T_mean_se(gi); ee = ee(order);
        ci = mod(g - 1, size(palette, 1)) + 1;
        mi = mod(g - 1, numel(markers)) + 1;
        errorbar(ax, vs, yy, ee, ee, markers{mi}, 'LineStyle', '-', ...
                 'Color', palette(ci, :), 'MarkerFaceColor', palette(ci, :), ...
                 'MarkerSize', 5, 'LineWidth', 1.3, 'CapSize', 4, ...
                 'DisplayName', group_labels{g});
    end

    % State-dependent subscript on the analytic value: T^{analytic}_r for the
    % (2,0,0) radial kinetic energy, T^{analytic}_theta for the (2,1,0) polar one.
    if strcmp(state_label, '2p0')
        T_an_label = '$T^{\mathrm{analytic}}_\theta$';
    else
        T_an_label = '$T^{\mathrm{analytic}}_r$';
    end

    set(ax, 'XScale', 'log');
    x_now = xlim(ax);
    plot(ax, x_now, [T_analytical, T_analytical], '-.', ...
         'Color', [0.60, 0.00, 0.00], 'LineWidth', 1.3, ...
         'DisplayName', T_an_label);
    xlim(ax, x_now);

    xlabel(ax, '$v_{\max}/c$', 'Interpreter', 'latex', 'FontSize', 14);
    if strcmp(state_label, '2p0')
        ylabel(ax, '$\langle T_\theta \rangle$ (J)', 'Interpreter', 'latex', ...
               'FontSize', 14);
        title(ax, '$(2,1,0)$ state', 'Interpreter', 'latex', 'FontSize', 14);
        ylim(ax, [2e-19, 9e-19]);
    else
        ylabel(ax, '$\langle T_r \rangle$ (J)', 'Interpreter', 'latex', ...
               'FontSize', 14);
        title(ax, '$(2,0,0)$ state', 'Interpreter', 'latex', 'FontSize', 14);
        ylim(ax, [2e-19, 9e-19]);
    end
    legend(ax, 'show', 'Interpreter', 'latex', 'Location', 'northwest', ...
           'Box', 'off', 'FontSize', 13);
    grid(ax, 'on');
    hold(ax, 'off');

    pdf_name = [state_label '_cutoff_scan_kinetic' name_suffix '_manuscript.pdf'];
    exportgraphics(fig, fullfile(export_dir, pdf_name), 'ContentType', 'vector');
    fprintf('%s written\n', fullfile(export_dir, pdf_name));
end

function [value, seed_id, dt_zs] = parse_vmax_from_filename(state_label, mat_name)
    % Parse v_max/c, seed, and dt (in zeptoseconds) from a scan .mat filename.
    % The seed is the last token; the dt tag is optional and an untagged run
    % uses the default dt = 10 zs (1e-20 s). Filenames look like:
    %   2s0_scan_vmax_1p0c_seed_1001.mat
    %   2s0_scan_vmax_5p0c_seed_112.mat
    %   2s0_scan_vmax_5p0c_dt_2zs_seed_10000.mat
    v_tokens = regexp(mat_name, [state_label '_scan.*_vmax_([0-9]+p[0-9]+c)'], ...
                      'tokens', 'once');
    if isempty(v_tokens)
        error('Could not parse v_max/c from filename: %s', mat_name);
    end
    value = str2double(strrep(strrep(v_tokens{1}, 'c', ''), 'p', '.'));

    seed_tokens = regexp(mat_name, '_seed_([0-9]+)', 'tokens', 'once');
    if isempty(seed_tokens)
        seed_id = NaN;
    else
        seed_id = str2double(seed_tokens{1});
    end

    dt_tokens = regexp(mat_name, '_dt_([0-9p]+)zs', 'tokens', 'once');
    if isempty(dt_tokens)
        dt_zs = 10;   % default dt = 1e-20 s = 10 zs
    else
        dt_zs = str2double(strrep(dt_tokens{1}, 'p', '.'));
    end
end

function plot_dt_grouped_errorbar(v_over_c, y_value, lower_err, upper_err, group_id, group_labels, marker)
    % One connected error-bar curve per group, x = v_max/c. Within a group the
    % points are sorted by v_max so the line joins them in order. Groups are by
    % dt (or by (dt, T) when mixing total times); see make_curve_groups.
    colors = line_colors(numel(group_labels));
    hold on;
    for g = 1:numel(group_labels)
        gi = (group_id == g);
        if ~any(gi), continue; end
        [vs, order] = sort(v_over_c(gi));
        yy = y_value(gi);  yy = yy(order);
        lo = lower_err(gi); lo = lo(order);
        hi = upper_err(gi); hi = hi(order);
        errorbar(vs, yy, lo, hi, marker, 'LineStyle', '-', ...
                 'Color', colors(g, :), 'MarkerFaceColor', colors(g, :), ...
                 'LineWidth', 1.3, 'DisplayName', group_labels{g});
    end
    hold off;
end

function plot_dt_grouped_line(v_over_c, y_value, group_id, group_labels, marker)
    % One connected marker-line per group (no error bars), x = v_max/c.
    colors = line_colors(numel(group_labels));
    hold on;
    for g = 1:numel(group_labels)
        gi = (group_id == g);
        if ~any(gi), continue; end
        [vs, order] = sort(v_over_c(gi));
        yy = y_value(gi);  yy = yy(order);
        plot(vs, yy, [marker '-'], 'Color', colors(g, :), ...
             'MarkerFaceColor', colors(g, :), 'LineWidth', 1.3, ...
             'DisplayName', group_labels{g});
    end
    hold off;
end

function plot_trajectory_points_panel(v_over_c, T_per_traj_cell, group_id, group_labels, T_analytical)
    % Scatter every per-trajectory time-averaged kinetic energy at x = v_max/c
    % (log axis). Points are grouped (colour) with a multiplicative per-group
    % x-offset and per-point jitter so the clouds separate, and the per-file
    % mean (diamond) and median (thick tick) are overlaid. Reuses the
    % per-trajectory vectors already loaded by the caller.
    ax = gca;
    hold(ax, 'on');
    set(ax, 'XScale', 'log');

    ng = numel(group_labels);
    colors = line_colors(ng);
    group_handles = gobjects(ng, 1);
    assigned = false(ng, 1);
    rng(1);   % reproducible jitter

    for gi = 1:ng
        col = colors(gi, :);
        % multiplicative x-offset per group so the clouds separate on the
        % log axis; the offset is symmetric about the true v_max/c.
        x_off = exp(0.08 * (gi - (ng + 1) / 2));
        idx = find(group_id == gi);
        for j = 1:numel(idx)
            k = idx(j);
            y = T_per_traj_cell{k};
            xc = v_over_c(k) * x_off;
            x = xc * exp(0.025 * randn(size(y)));      % per-point jitter
            h = scatter(ax, x, y, 10, col, 'filled', ...
                        'MarkerFaceAlpha', 0.30, 'MarkerEdgeColor', 'none', ...
                        'HandleVisibility', 'off');
            if ~assigned(gi)
                group_handles(gi) = h;
                assigned(gi) = true;
                set(h, 'HandleVisibility', 'on', 'DisplayName', group_labels{gi});
            end
            % per-file mean (diamond) and median (thick tick)
            plot(ax, xc, mean(y), 'd', 'MarkerSize', 7, ...
                 'MarkerFaceColor', col, 'MarkerEdgeColor', 'k', ...
                 'LineWidth', 0.8, 'HandleVisibility', 'off');
            plot(ax, xc * [1/1.06, 1.06], median(y) * [1, 1], '-', ...
                 'Color', max(col - 0.15, 0), 'LineWidth', 2.2, ...
                 'HandleVisibility', 'off');
        end
    end

    xlabel(ax, '$v_{\max}/c$', 'Interpreter', 'latex');
    ylabel(ax, '$\left< T \right>_j$ per trajectory (J)', 'Interpreter', 'latex');
    grid(ax, 'on');
    draw_horizontal_line(T_analytical, '-.', [0.60, 0.00, 0.00]);
    legend(ax, 'show', 'Interpreter', 'latex', 'Location', 'best', 'Box', 'off');
    hold(ax, 'off');
end

function require_field(S, field_name, mat_name)
    if ~isfield(S, field_name)
        error('Required field "%s" is missing from %s.', field_name, mat_name);
    end
end

function T_time = optional_time_series(S, field_name)
    if isfield(S, field_name)
        T_time = S.(field_name);
        if isempty(T_time) || all(T_time(:) == 0)
            T_time = [];
        end
    else
        T_time = [];
    end
end

function T_analytical = kinetic_energy_reference_value(state_label)
    if strcmp(state_label, '2p0')
        T_analytical = 3.6330e-19;
    else
        T_analytical = 5.4497e-19;
    end
end

function s = dt_seconds_latex(dt_zs)
    % Format a dt given in zeptoseconds on a common 10^{-21} s basis so the
    % values are directly comparable in the legend, e.g.
    % 100 zs -> '100 \times 10^{-21}', 10 zs -> '10 \times 10^{-21}',
    % 5 zs -> '5 \times 10^{-21}', 0.5 zs -> '0.5 \times 10^{-21}'.
    s = sprintf('%g \\times 10^{-21}', dt_zs);
end

function s = seconds_sci_latex(t)
    % Format a time in seconds as a compact LaTeX mantissa x 10^exp string,
    % e.g. 2e-11 -> '2\times10^{-11}', 1e-11 -> '10^{-11}'.
    if ~isfinite(t) || t <= 0
        s = '?'; return;
    end
    e = floor(log10(t) + 1e-9);
    m = round(t / 10^e * 1e6) / 1e6;
    if m >= 10 - 1e-6
        m = m / 10; e = e + 1;
    end
    if abs(m - 1) < 1e-9
        s = sprintf('10^{%d}', e);
    else
        s = sprintf('%g \\times 10^{%d}', m, e);
    end
end

function [group_id, group_labels] = make_curve_groups(dt_zs, total_time, mix_time)
    % Assign each scan point to a curve group and build its legend label.
    % Default: group by integration step dt only (legend shows dt). With
    % mix_time on: group by (dt, total time T) so runs of different total time
    % stay on separate, correctly labelled curves (legend shows dt and T).
    dt_zs = dt_zs(:);
    if mix_time
        % Round T to 1e-13 s to fold floating-point-identical grids together.
        t_round = round(total_time(:) / 1e-13);
        [keys, ~, group_id] = unique([dt_zs, t_round], 'rows');  % sorted by dt then T
        group_labels = cell(1, size(keys, 1));
        for g = 1:size(keys, 1)
            group_labels{g} = sprintf('$\\Delta t = %s$ s, $T = %s$ s', ...
                dt_seconds_latex(keys(g, 1)), seconds_sci_latex(keys(g, 2) * 1e-13));
        end
    else
        [keys, ~, group_id] = unique(dt_zs);   % sorted ascending
        group_labels = cell(1, numel(keys));
        for g = 1:numel(keys)
            group_labels{g} = sprintf('$\\Delta t = %s$ s', dt_seconds_latex(keys(g)));
        end
    end
end

function line_style = line_style_for_seed(seed_value)
    if isnan(seed_value)
        line_style = '-';
    else
        line_style = '--';
    end
end

function colors = line_colors(n)
    base = [ ...
        0.00, 0.25, 0.55; ...
        0.10, 0.45, 0.25; ...
        0.55, 0.25, 0.00; ...
        0.45, 0.10, 0.55; ...
        0.50, 0.30, 0.05; ...
        0.10, 0.55, 0.60];
    if n <= size(base, 1)
        colors = base(1:n, :);
    else
        colors = lines(n);
    end
end

function draw_horizontal_line(y_value, line_style, color)
    x_limits = xlim;
    plot(x_limits, [y_value, y_value], line_style, 'Color', color, ...
         'LineWidth', 1.2, 'HandleVisibility', 'off');
    xlim(x_limits);
end

function se = bootstrap_mean_se(values, n_boot, seed)
    % Standard error of the sample mean estimated by nonparametric
    % bootstrap. The classical SEM = std(x)/sqrt(n) underestimates the
    % seed-to-seed uncertainty when the underlying distribution has a
    % heavy (Lomax) tail, because the empirical variance is bounded by
    % the relativistic cutoff rather than by the true (infinite)
    % population variance. The bootstrap captures the rare-event
    % contribution honestly.
    n = numel(values);
    if n < 2 || n_boot < 2
        se = NaN; return;
    end
    rng_state = rng();
    cleanup = onCleanup(@() rng(rng_state));   %#ok<NASGU>
    rng(seed);
    boot_idx = randi(n, n, n_boot);
    means = mean(values(boot_idx), 1);
    se = std(means);
end

function se = bootstrap_median_se(values, n_boot, seed)
    % Standard error of the sample median estimated by nonparametric
    % bootstrap. The IQR (Q3 - Q1) describes the spread of the empirical
    % within-sample distribution and converges to a finite population
    % IQR as M -> Inf; it is NOT the uncertainty of the median estimator.
    n = numel(values);
    if n < 2 || n_boot < 2
        se = NaN; return;
    end
    rng_state = rng();
    cleanup = onCleanup(@() rng(rng_state));   %#ok<NASGU>
    rng(seed);
    boot_idx = randi(n, n, n_boot);
    medians = median(values(boot_idx), 1);
    se = std(medians);
end

function plot_light_whiskers(v_over_c, y_value, lower_err, upper_err)
    % Thin light-gray whiskers used as a descriptive overlay (classical
    % SEM on panel (a), IQR on panel (b)). They are not the primary
    % error bars and are intentionally drawn without markers so they do
    % not compete with the bootstrap-SE error bars on the same axes. The
    % grouping does not affect their appearance, so all points are drawn
    % with a single gray errorbar call.
    light_gray = [0.55, 0.55, 0.55];
    hold on;
    errorbar(v_over_c, y_value, lower_err, upper_err, ...
             'LineStyle', 'none', 'Color', light_gray, ...
             'Marker', 'none', 'LineWidth', 0.6, 'CapSize', 4, ...
             'HandleVisibility', 'off');
    hold off;
end
