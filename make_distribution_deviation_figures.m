function make_distribution_deviation_figures(state)
% MAKE_DISTRIBUTION_DEVIATION_FIGURES  Build Fig. 8 deviation panel for one state.
%
% Usage:
%   make_distribution_deviation_figures("2p0")
%   make_distribution_deviation_figures("2s0")

    addpath('functions');
    if nargin < 1 || isempty(state)
        error('Usage: make_distribution_deviation_figures("2p0") or "2s0"');
    end

    cfg = deviation_config(state);
    require_file(cfg.mat_path);
    S = load(cfg.mat_path);

    fig = manuscript_figure(20);
    plot_distributions_deviation(S.time_arr, S.dev_arr_R, S.dev_arr_theta, S.dev_arr_phi);
    manuscript_bump_fonts(18);
    exportgraphics(fig, fullfile('figures', cfg.pdf_name), 'ContentType', 'vector');
    fprintf('%s -> t=%.1f fs\n', cfg.pdf_name, S.time_arr(1, end) / 1e-15);
    close(fig);
end

function cfg = deviation_config(state)
    key = lower(strtrim(char(state)));
    key = strrep(key, '-', '_');
    switch key
        case '2p0'
            cfg.mat_path = fullfile('data', '2p0_10ps', '2p0_10ps.mat');
            cfg.pdf_name = '2p0_10ps-distribution_deviation.pdf';
        case '2s0'
            cfg.mat_path = fullfile('data', '2s0_10ps', '2s0_10ps.mat');
            cfg.pdf_name = '2s0_10ps-distribution_deviation.pdf';
        otherwise
            error('Unknown state "%s". Use "2p0" or "2s0".', state);
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
