function make_distribution_figures(state)
% MAKE_DISTRIBUTION_FIGURES  Build manuscript distribution panels for one state.
%
% Usage:
%   make_distribution_figures("1s0")
%   make_distribution_figures("2p0")
%   make_distribution_figures("2p_m1")

    addpath('functions');
    if nargin < 1 || isempty(state)
        error('Usage: make_distribution_figures("1s0"), "2p0", or "2p_m1"');
    end

    cfg = distribution_config(state);
    require_file(cfg.mat_path);
    S = load(cfg.mat_path);
    [P_r, P_theta, P_phi] = manuscript_analytic_distributions(S);

    fig = manuscript_figure(20);
    for k = 1:numel(cfg.panels)
        panel = cfg.panels(k);
        clf(fig);
        switch panel.kind
            case 'radial_traj'
                plot_radial_distribution(S.hist_bins_r, S.hist_counts_r_traj(panel.index, :), ...
                                         S.r_values, P_r);
            case 'polar_traj'
                plot_polar_distribution(S.hist_bins_theta, S.hist_counts_theta_traj(panel.index, :), ...
                                        S.theta_values, P_theta);
            case 'azimuthal_traj'
                plot_azimuthal_distribution(S.hist_bins_phi, S.hist_counts_phi_traj(panel.index, :), ...
                                            S.phi_values, P_phi);
            case 'radial_full'
                plot_radial_distribution(S.hist_bins_r, S.hist_counts_r, S.r_values, P_r);
            case 'polar_full'
                plot_polar_distribution(S.hist_bins_theta, S.hist_counts_theta, S.theta_values, P_theta);
            case 'azimuthal_full'
                plot_azimuthal_distribution(S.hist_bins_phi, S.hist_counts_phi, S.phi_values, P_phi);
            otherwise
                error('Unknown distribution panel kind "%s".', panel.kind);
        end
        manuscript_bump_fonts(20);
        add_panel_time_label(fig, panel_time_label(S, panel.index));
        exportgraphics(fig, fullfile('figures', panel.pdf_name), 'ContentType', 'vector');
        fprintf('%s -> %s\n', panel.pdf_name, panel.time_desc);
    end
    close(fig);
end

function cfg = distribution_config(state)
    key = normalize_state(state);
    switch key
        case '1s0'
            cfg.mat_path = fullfile('data', '1s0_1ps', '1s0_1ps.mat');
            cfg.panels = [
                panel('radial_traj', 500, '1s0_1ps-radial_hist_500fs.pdf', 't=500 fs')
                panel('polar_traj', 500, '1s0_1ps-polar_hist_500fs.pdf', 't=500 fs')
                panel('azimuthal_traj', 500, '1s0_1ps-azimuthal_hist_500fs.pdf', 't=500 fs')
            ];

        case '2p0'
            cfg.mat_path = fullfile('data', '2p0_10ps', '2p0_10ps.mat');
            cfg.panels = [
                panel('radial_traj', 500, '2p0_10ps-radial_hist_5ps.pdf', 't=5 ps')
                panel('polar_traj', 500, '2p0_10ps-polar_hist_5ps.pdf', 't=5 ps')
                panel('azimuthal_traj', 500, '2p0_10ps-azimuthal_hist_5ps.pdf', 't=5 ps')
            ];

        case '2p_m1'
            cfg.mat_path = fullfile('data', '2p_m1_1ps', '2p_m1_1ps.mat');
            cfg.panels = [
                panel('radial_full', NaN, '2p_m1_1ps-radial_hist.pdf', 'full run')
                panel('polar_full', NaN, '2p_m1_1ps-polar_hist.pdf', 'full run')
                panel('azimuthal_full', NaN, '2p_m1_1ps-azimuthal_hist.pdf', 'full run')
            ];

        otherwise
            error('Unknown state "%s". Use "1s0", "2p0", or "2p_m1".', state);
    end
end

function p = panel(kind, index, pdf_name, time_desc)
    p = struct('kind', kind, 'index', index, 'pdf_name', pdf_name, 'time_desc', time_desc);
end

function key = normalize_state(state)
    key = lower(strtrim(char(state)));
    key = strrep(key, '-', '_');
    if any(strcmp(key, {'2p1', '2p_m+1', '2p_p1'}))
        key = '2p_m1';
    end
end

function require_file(path)
    if ~exist(path, 'file')
        error('Missing %s. Run the corresponding simulation first.', path);
    end
    if ~exist('figures', 'dir')
        mkdir('figures');
    end
end

function label = panel_time_label(S, index)
    if isnan(index)
        seconds = S.time_arr(end);
    else
        seconds = S.time_arr(index);
    end
    label = format_time_label(seconds);
end

function label = format_time_label(seconds)
    if ~isfinite(seconds)
        label = 't = n/a';
    elseif seconds < 1e-12
        label = sprintf('t = %.0f fs', seconds / 1e-15);
    else
        label = sprintf('t = %.0f ps', seconds / 1e-12);
    end
end

function add_panel_time_label(fig, label)
    annotation(fig, 'textbox', [0.15 0.925 0.34 0.06], ...
               'String', label, 'EdgeColor', 'none', ...
               'Color', [0.10 0.10 0.10], ...
               'FontSize', 24, 'FontWeight', 'normal', ...
               'Interpreter', 'tex', 'VerticalAlignment', 'bottom');
end
