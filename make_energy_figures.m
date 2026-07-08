function make_energy_figures(state)
% MAKE_ENERGY_FIGURES  Build manuscript energy panel for one state.
%
% Usage:
%   make_energy_figures("2p0")
%   make_energy_figures("2s0")

    addpath('functions');
    if nargin < 1 || isempty(state)
        error('Usage: make_energy_figures("2p0") or "2s0"');
    end

    cfg = energy_config(state);
    require_file(cfg.mat_path);
    S = load(cfg.mat_path);

    if isfield(cfg, 'index_fraction')
        idx = 1:size(S.time_arr, 2) / cfg.index_fraction;
    else
        idx = 1:size(S.time_arr, 2);
    end

    fig = manuscript_figure(20);
    plot_energies_with_bands(S.n, S.l, S.m, S.time_arr(idx), ...
        S.KE_radial_traj_M(:, idx), S.KE_theta_traj_M(:, idx), S.KE_phi_traj_M(:, idx), ...
        S.V_traj_M(:, idx), S.E_traj_M(:, idx));
    manuscript_bump_fonts(18);
    exportgraphics(fig, fullfile('figures', cfg.pdf_name), 'ContentType', 'vector');
    fprintf('%s -> t=%.1f fs\n', cfg.pdf_name, S.time_arr(1, idx(end)) / 1e-15);
    close(fig);
end

function cfg = energy_config(state)
    key = lower(strtrim(char(state)));
    key = strrep(key, '-', '_');
    switch key
        case '2p0'
            cfg.mat_path = fullfile('data', '2p0_10ps', '2p0_10ps.mat');
            cfg.pdf_name = '2p0_1ps-energies.pdf';
            cfg.index_fraction = 10;
        case '2s0'
            cfg.mat_path = fullfile('data', '2s0_10ps', '2s0_10ps.mat');
            cfg.pdf_name = '2s0_10ps-energies.pdf';
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
