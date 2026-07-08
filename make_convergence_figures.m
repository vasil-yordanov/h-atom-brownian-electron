function make_convergence_figures(state)
% MAKE_CONVERGENCE_FIGURES  Build manuscript convergence panels for one state.
%
% Usage:
%   make_convergence_figures("2p0")
%   make_convergence_figures("2s0")

    addpath('functions');
    if nargin < 1 || isempty(state)
        error('Usage: make_convergence_figures("2p0") or "2s0"');
    end

    cfg = convergence_config(state);
    require_file(cfg.mat_path);
    S = load(cfg.mat_path);
    [P_r, P_theta, ~] = manuscript_analytic_distributions(S);

    fig = manuscript_figure(20);
    for k = 1:numel(cfg.panels)
        panel = cfg.panels(k);
        clf(fig);
        switch panel.kind
            case 'polar'
                plot_polar_distribution(S.hist_bins_theta, S.hist_counts_theta_traj(panel.index, :), ...
                                        S.theta_values, P_theta);
            case 'radial'
                plot_radial_distribution(S.hist_bins_r, S.hist_counts_r_traj(panel.index, :), ...
                                         S.r_values, P_r);
            otherwise
                error('Unknown convergence panel kind "%s".', panel.kind);
        end
        manuscript_bump_fonts(20);
        add_panel_time_label(fig, panel_time_label(S, panel.index));
        exportgraphics(fig, fullfile('figures', panel.pdf_name), 'ContentType', 'vector');
        fprintf('%s -> %s\n', panel.pdf_name, panel.time_desc);
    end
    close(fig);
end

function cfg = convergence_config(state)
    key = lower(strtrim(char(state)));
    key = strrep(key, '-', '_');
    switch key
        case '2p0'
            cfg.mat_path = fullfile('data', '2p0_10ps', '2p0_10ps.mat');
            cfg.panels = [
                panel('polar', 10, '2p0_10ps-polar_hist_100fs.pdf', 't=100 fs')
                panel('polar', 100, '2p0_10ps-polar_hist_1ps.pdf', 't=1 ps')
                panel('polar', 1000, '2p0_10ps-polar_hist_10ps.pdf', 't=10 ps')
            ];

        case '2s0'
            cfg.mat_path = fullfile('data', '2s0_10ps', '2s0_10ps.mat');
            cfg.panels = [
                panel('radial', 60, '2s0_10ps-radial_hist_300fs.pdf', 't=300 fs')
                panel('radial', 68, '2s0_10ps-radial_hist_340fs.pdf', 't=340 fs')
                panel('radial', 100, '2s0_10ps-radial_hist_500fs.pdf', 't=500 fs')
            ];

        otherwise
            error('Unknown state "%s". Use "2p0" or "2s0".', state);
    end
end

function p = panel(kind, index, pdf_name, time_desc)
    p = struct('kind', kind, 'index', index, 'pdf_name', pdf_name, 'time_desc', time_desc);
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
    seconds = S.time_arr(index);
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
