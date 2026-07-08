function varargout = make_lz_running_estimator(state_name, figure_dir)
% MAKE_LZ_RUNNING_ESTIMATOR  Running Lz estimator from stored trajectories.
%
% Usage:
%   make_lz_running_estimator
%   make_lz_running_estimator("2p_m1")
%   make_lz_running_estimator("2p0")
%   make_lz_running_estimator("2p_mn1")
%
% The script prefers compact online Lz estimator runs (*_lz_1ps). If those are
% absent, it falls back to post-processing stored r_traj/theta_traj/phi_traj
% arrays. The post-processed estimator uses only kinematic trajectory data:
%
%   Lz(t) = cumsum(m_e*r^2*sin(theta)^2*dphi) / t
%
% with unwrap(phi) before diff(phi), and start-of-interval weights.
%
% Output:
%   figures/lz_running_estimator.pdf

    addpath('functions');
    if nargin < 1 || isempty(state_name)
        state_name = "all";
    end
    if nargin < 2 || isempty(figure_dir)
        figure_dir = 'figures';
    end
    if ~exist(figure_dir, 'dir')
        mkdir(figure_dir);
    end

    [hbar, m_e] = constants();
    configs = select_configs(string(state_name));
    results = repmat(empty_result(), 0, 1);

    for k = 1:numel(configs)
        cfg = resolve_config(configs(k));
        fprintf('\nState %s: %s\n', cfg.state_label, cfg.mat_path);
        result = process_one_state(cfg, hbar, m_e);
        results(end + 1, 1) = result; %#ok<AGROW>
    end

    print_summary(results);
    make_lz_figure(results, figure_dir);
    if nargout > 0
        varargout{1} = results;
    end
end

function configs = select_configs(state_name)
    all_configs = lz_configs();
    if any(strcmpi(state_name, ["all", "2p", "2p_all"]))
        configs = all_configs;
        return;
    end

    aliases = strings(size(all_configs));
    for k = 1:numel(all_configs)
        aliases(k) = all_configs(k).key;
    end
    hit = false(size(all_configs));
    for k = 1:numel(all_configs)
        hit(k) = any(strcmpi(state_name, all_configs(k).aliases));
    end
    if ~any(hit)
        error('Unknown state "%s". Use "all", "2p_m1", "2p0", or "2p_mn1".', state_name);
    end
    configs = all_configs(hit);
end

function configs = lz_configs()
    configs(1) = struct( ...
        'key', "2p_p1", ...
        'aliases', ["2p_p1", "2p_plus1", "2p_mplus1", "2p_m1", "2p1", "+1", "m=+1"], ...
        'state_label', "(2,1,+1)", ...
        'legend_label', "$m=+1$", ...
        'candidates', ["2p_m1_lz_1ps", "2p_m1_movie", "2p_m1_1ps"]);
    configs(2) = struct( ...
        'key', "2p_0", ...
        'aliases', ["2p_0", "2p0", "2p_m0", "0", "m=0"], ...
        'state_label', "(2,1,0)", ...
        'legend_label', "$m=\phantom{+}0$", ...
        'candidates', ["2p_m0_lz_1ps", "2p_m0_movie", "2p0_movie", "2p_m0_1ps"]);
    configs(3) = struct( ...
        'key', "2p_mn1", ...
        'aliases', ["2p_mn1", "2p_minus1", "2p_m-1", "-1", "m=-1"], ...
        'state_label', "(2,1,-1)", ...
        'legend_label', "$m=-1$", ...
        'candidates', ["2p_mn1_lz_1ps", "2p_mn1_movie", "2p_mn1_1ps"]);
end

function cfg = resolve_config(cfg)
    for k = 1:numel(cfg.candidates)
        run_name = char(cfg.candidates(k));
        mat_path = fullfile('data', run_name, [run_name '.mat']);
        if exist(mat_path, 'file')
            cfg.run_name = string(run_name);
            cfg.mat_path = mat_path;
            cfg.uses_lz_run = contains(cfg.run_name, "_lz_");
            cfg.uses_movie_run = endsWith(cfg.run_name, "_movie");
            if ~cfg.uses_lz_run && ~cfg.uses_movie_run
                warning(['Using %s because no *_lz_1ps or *_movie dataset was found for %s. ', ...
                         'This may contain only the first stored part of the run.'], ...
                        cfg.run_name, cfg.state_label);
            end
            return;
        end
    end
    error(['No .mat file found for %s. Expected one of: %s. ', ...
           'For the 1 ps Lz diagnostic, run the matching *_lz_1ps preset first.'], ...
          cfg.state_label, strjoin(cfg.candidates, ', '));
end

function result = process_one_state(cfg, hbar, m_e)
    saved_names = string(who('-file', cfg.mat_path));
    if any(saved_names == "Lz_over_hbar_traj")
        result = process_online_lz_state(cfg);
        return;
    end

    meta = load(cfg.mat_path, 'dt', 'n_steps');
    mf = matfile(cfg.mat_path);
    n_stored = stored_count(mf, 'r_traj');
    if n_stored < 2
        error('%s has fewer than two stored trajectory points.', cfg.mat_path);
    end

    dts = meta.dt;
    decimation_factor = dts / meta.dt;
    if decimation_factor >= 1e3
        error(['Stored decimation factor K=%.4g is too coarse for %s. ', ...
               'Add online accumulation in h_atom.m and rerun.'], ...
              decimation_factor, cfg.state_label);
    elseif decimation_factor > 100
        warning('Stored decimation factor K=%.4g is above the preferred K<=100 range.', ...
                decimation_factor);
    end

    n_inc = n_stored - 1;
    n_out = min(5000, n_inc);
    out_idx = unique(round(linspace(1, n_inc, n_out)));
    mid_idx = max(1, round(n_inc / 2));
    robust_inc = floor((n_stored - 1) / 2);
    robust_full_idx = 2 * robust_inc;

    if n_stored <= 3e7
        traj = load(cfg.mat_path, 'r_traj', 'theta_traj', 'phi_traj');
        r_all = traj.r_traj(:);
        th_all = traj.theta_traj(:);
        ph_all = traj.phi_traj(:);
        full = memory_estimator(r_all, th_all, ph_all, dts, 1, out_idx, ...
                                unique([mid_idx, n_inc, robust_full_idx]), hbar, m_e);
        robust = memory_estimator(r_all, th_all, ph_all, dts, 2, [], robust_inc, hbar, m_e);
    else
        full = stream_estimator(mf, n_stored, dts, 1, out_idx, ...
                                unique([mid_idx, n_inc, robust_full_idx]), hbar, m_e);
        robust = stream_estimator(mf, n_stored, dts, 2, [], robust_inc, hbar, m_e);
    end

    result = empty_result();
    result.key = cfg.key;
    result.state_label = cfg.state_label;
    result.legend_label = cfg.legend_label;
    result.run_name = cfg.run_name;
    result.mat_path = string(cfg.mat_path);
    result.uses_movie_run = cfg.uses_movie_run;
    result.dt = meta.dt;
    result.dts = dts;
    result.decimation_factor = decimation_factor;
    result.n_stored = n_stored;
    result.t_final = n_inc * dts;
    result.lz_mid = full.values(mid_idx);
    result.lz_final = full.values(n_inc);
    result.robust_t_final = robust_inc * 2 * dts;
    result.robust_lz_final = robust.values(robust_inc);
    result.robust_full_lz = full.values(robust_full_idx);
    result.robust_diff = result.robust_lz_final - result.robust_full_lz;
    result.max_abs_dphi = full.max_abs_dphi;
    result.t_plot = full.t_out(:);
    result.lz_plot = full.lz_out(:);

    fprintf('  stored samples        : %d\n', n_stored);
    fprintf('  dt, dts, K            : %.3g s, %.3g s, %.4g\n', ...
            result.dt, result.dts, result.decimation_factor);
    fprintf('  stored trajectory T   : %.6g ps\n', result.t_final / 1e-12);
    fprintf('  max |dphi| stored     : %.4g rad\n', result.max_abs_dphi);
end

function result = process_online_lz_state(cfg)
    S = load(cfg.mat_path, 'time_arr', 'Lz_time_arr', 'Lz_over_hbar_traj', ...
             'dt', 'frame_step');
    if isfield(S, 'Lz_time_arr')
        t = S.Lz_time_arr(:);
    else
        t = S.time_arr(:);
    end
    lz = S.Lz_over_hbar_traj(:);
    valid = isfinite(lz);
    t = t(valid);
    lz = lz(valid);
    if isempty(t)
        error('%s has Lz_over_hbar_traj but no finite samples.', cfg.mat_path);
    end

    mid_idx = max(1, round(numel(t) / 2));
    result = empty_result();
    result.key = cfg.key;
    result.state_label = cfg.state_label;
    result.legend_label = cfg.legend_label;
    result.run_name = cfg.run_name;
    result.mat_path = string(cfg.mat_path);
    result.method = "online";
    result.uses_lz_run = true;
    result.uses_movie_run = false;
    result.dt = S.dt;
    result.dts = S.frame_step * S.dt;
    result.decimation_factor = S.frame_step;
    result.n_stored = numel(t);
    result.t_final = t(end);
    result.lz_mid = lz(mid_idx);
    result.lz_final = lz(end);
    result.robust_t_final = NaN;
    result.robust_lz_final = NaN;
    result.robust_full_lz = NaN;
    result.robust_diff = NaN;
    result.max_abs_dphi = NaN;
    result.t_plot = t;
    result.lz_plot = lz;

    fprintf('  method                : online accumulation\n');
    fprintf('  saved samples         : %d\n', numel(t));
    fprintf('  dt, output interval   : %.3g s, %.3g s\n', result.dt, result.dts);
    fprintf('  saved trajectory T    : %.6g ps\n', result.t_final / 1e-12);
end

function out = stream_estimator(mf, n_stored, dts, stride, out_idx, value_idx, hbar, m_e)
    n_inc = floor((n_stored - 1) / stride);
    chunk_inc = 5e6;
    running_sum = 0;
    max_abs_dphi = 0;

    out.t_out = zeros(numel(out_idx), 1);
    out.lz_out = zeros(numel(out_idx), 1);
    out.values = containers.Map('KeyType', 'double', 'ValueType', 'double');
    out.max_abs_dphi = 0;

    out_pos = 1;
    value_idx = value_idx(:)';
    inc_start = 1;
    while inc_start <= n_inc
        inc_end = min(n_inc, inc_start + chunk_inc - 1);
        start_samples = 1 + (inc_start - 1) * stride;
        end_samples = 1 + (inc_end - 1) * stride;
        sample_idx = start_samples:stride:(end_samples + stride);

        r = read_traj(mf, 'r_traj', sample_idx(1:end - 1));
        th = read_traj(mf, 'theta_traj', sample_idx(1:end - 1));
        ph = read_traj(mf, 'phi_traj', sample_idx);

        phiu = unwrap(ph(:));
        dphi = diff(phiu);
        max_abs_dphi = max(max_abs_dphi, max(abs(dphi)));

        weights = m_e .* r(:).^2 .* sin(th(:)).^2;
        partial = running_sum + cumsum(weights .* dphi);
        while out_pos <= numel(out_idx) && out_idx(out_pos) <= inc_end
            if out_idx(out_pos) >= inc_start
                loc = out_idx(out_pos) - inc_start + 1;
                out.t_out(out_pos) = out_idx(out_pos) * stride * dts;
                out.lz_out(out_pos) = partial(loc) ./ (out_idx(out_pos) * stride * dts) ./ hbar;
            end
            out_pos = out_pos + 1;
        end

        wanted = value_idx(value_idx >= inc_start & value_idx <= inc_end);
        for w = wanted
            loc = w - inc_start + 1;
            out.values(w) = partial(loc) ./ (w * stride * dts) ./ hbar;
        end

        running_sum = partial(end);
        inc_start = inc_end + 1;
    end

    out.max_abs_dphi = max_abs_dphi;
end

function out = memory_estimator(r, th, ph, dts, stride, out_idx, value_idx, hbar, m_e)
    idx = 1:stride:numel(ph);
    n_inc = numel(idx) - 1;
    phiu = unwrap(ph(idx));
    dphi = diff(phiu);
    weights = m_e .* r(idx(1:end - 1)).^2 .* sin(th(idx(1:end - 1))).^2;
    partial = cumsum(weights(:) .* dphi(:));

    out.t_out = zeros(numel(out_idx), 1);
    out.lz_out = zeros(numel(out_idx), 1);
    out.values = containers.Map('KeyType', 'double', 'ValueType', 'double');
    out.max_abs_dphi = max(abs(dphi));

    for k = 1:numel(out_idx)
        ii = out_idx(k);
        out.t_out(k) = ii * stride * dts;
        out.lz_out(k) = partial(ii) ./ (ii * stride * dts) ./ hbar;
    end

    value_idx = value_idx(:)';
    value_idx = value_idx(value_idx >= 1 & value_idx <= n_inc);
    for ii = value_idx
        out.values(ii) = partial(ii) ./ (ii * stride * dts) ./ hbar;
    end
end

function n = stored_count(mf, var_name)
    sz = size(mf, var_name);
    n = max(sz);
end

function v = read_traj(mf, var_name, idx)
    sz = size(mf, var_name);
    if sz(1) == 1
        switch var_name
            case 'r_traj'
                v = mf.r_traj(1, idx);
            case 'theta_traj'
                v = mf.theta_traj(1, idx);
            case 'phi_traj'
                v = mf.phi_traj(1, idx);
            otherwise
                error('Unknown trajectory variable "%s".', var_name);
        end
    else
        switch var_name
            case 'r_traj'
                v = mf.r_traj(idx, 1);
            case 'theta_traj'
                v = mf.theta_traj(idx, 1);
            case 'phi_traj'
                v = mf.phi_traj(idx, 1);
            otherwise
                error('Unknown trajectory variable "%s".', var_name);
        end
    end
end

function print_summary(results)
    fprintf('\nRunning Lz estimator summary\n');
    fprintf('%-6s %-18s %-11s %8s %10s %12s %12s\n', ...
            'm', 'run', 'method', 'K', 'T(ps)', 'T/2', 'T');
    for k = 1:numel(results)
        fprintf('%-6s %-18s %-11s %8.4g %10.4g %12.6f %12.6f\n', ...
                results(k).legend_label, results(k).run_name, results(k).method, ...
                results(k).decimation_factor, ...
                results(k).t_final / 1e-12, results(k).lz_mid, ...
                results(k).lz_final);
        if isfinite(results(k).robust_diff)
            fprintf('       every-2 check at T: %+g\n', results(k).robust_diff);
        end
    end
end

function make_lz_figure(results, figure_dir)
    if isempty(results)
        return;
    end
    fig = manuscript_figure(20);
    set(fig, 'Units', 'centimeters', 'Position', [2, 2, 20, 14]);
    ax = axes('Parent', fig);
    hold(ax, 'on');
    colors = [0.8500, 0.3300, 0.1000;
              0.3500, 0.3500, 0.3500;
              0.0000, 0.4500, 0.7400];
    for k = 1:numel(results)
        plot(ax, results(k).t_plot, results(k).lz_plot, ...
             'LineWidth', 1.8, 'Color', colors(k, :), ...
             'DisplayName', char(results(k).legend_label));
    end
    yline(ax, 1, '--k', 'HandleVisibility', 'off');
    yline(ax, 0, '--k', 'HandleVisibility', 'off');
    yline(ax, -1, '--k', 'HandleVisibility', 'off');
    % Analytical +/-1 standard-deviation envelope for the m = +/-1 states:
    % std/hbar = sqrt(24*m_e*a0^2/(hbar*t)), from <r^2 sin^2 theta> = 24 a0^2
    % for the 2p, m = +/-1 Born density (manuscript Eq. for std[Lz_bar(t)]).
    [hbar_env, m_e_env] = constants();
    a0 = 5.29177210903e-11;
    t_env = linspace(results(1).t_plot(2), results(1).t_plot(end), 500);
    std_env = sqrt(24 * m_e_env * a0^2 ./ (hbar_env * t_env));
    for base = [1, -1]
        plot(ax, t_env, base + std_env, ':', 'Color', [0.35, 0.35, 0.35], ...
             'LineWidth', 0.9, 'HandleVisibility', 'off');
        plot(ax, t_env, base - std_env, ':', 'Color', [0.35, 0.35, 0.35], ...
             'LineWidth', 0.9, 'HandleVisibility', 'off');
    end
    ylim(ax, [-1.5, 1.5]);
    xlabel(ax, 'Time (s)');
    ylabel(ax, '$\bar{L}_z / \hbar$', 'Interpreter', 'latex');
    set(ax, 'FontSize', 18, 'Box', 'on', 'LineWidth', 1.0);
    ax.XAxis.Exponent = -12;
    ax.XLabel.FontSize = 20;
    ax.YLabel.FontSize = 20;
    lg = legend(ax, 'Interpreter', 'latex', 'Box', 'off', 'FontSize', 16);
    lg.Units = 'normalized';
    lg.Position(1:2) = [0.70, 0.58];
    grid(ax, 'on');
    hold(ax, 'off');

    pdf_path = fullfile(figure_dir, 'lz_running_estimator.pdf');
    exportgraphics(fig, pdf_path, 'ContentType', 'vector');
    close(fig);
    fprintf('Figure: %s\n', pdf_path);
end

function result = empty_result()
    result = struct( ...
        'key', "", ...
        'state_label', "", ...
        'legend_label', "", ...
        'run_name', "", ...
        'mat_path', "", ...
        'method', "postprocess", ...
        'uses_lz_run', false, ...
        'uses_movie_run', false, ...
        'dt', NaN, ...
        'dts', NaN, ...
        'decimation_factor', NaN, ...
        'n_stored', NaN, ...
        't_final', NaN, ...
        'lz_mid', NaN, ...
        'lz_final', NaN, ...
        'robust_t_final', NaN, ...
        'robust_lz_final', NaN, ...
        'robust_full_lz', NaN, ...
        'robust_diff', NaN, ...
        'max_abs_dphi', NaN, ...
        't_plot', [], ...
        'lz_plot', []);
end
