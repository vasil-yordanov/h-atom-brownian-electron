clearvars -except params_set_name extended;

% extended (default false): also generate the diagnostic (2,1,m) figures
% (per-state energies and xy clouds for the m = 0, +1, -1 states). These are not
% used in the manuscript .tex or the README, so they are off by default. Enable
% with:  extended = true; h_atom_post_processing
if ~exist('extended', 'var') || isempty(extended)
    extended = false;
end

export_base_path = 'figures';
load_base_path = 'data';

status = mkdir(export_base_path);

if status
    disp(['Directory "', export_base_path, '" created successfully.']);
else
    disp(['Failed to create directory "', export_base_path, '".']);
    exit;
end

if exist('params_set_name', 'var') && is_scan_parameter_set(params_set_name)
    maybe_analyse_cutoff_scan(params_set_name, export_base_path);
    return;
end

fig = figure;
defaultFontSize = 20;

% Set default font properties for the figure
set(fig, 'DefaultAxesFontSize', defaultFontSize, ...
         'DefaultTextFontSize', defaultFontSize, ...
         'DefaultLegendFontSize', defaultFontSize);

%%%%%%%%%%%%%%%%%%%%%

load(fullfile(load_base_path, '1s0_1fs/1s0_1fs.mat')); % fig_1
recompute_analytic_distributions;

% Manuscript Fig. 1: the 1s electron's partial trajectory after 200 as and 500 as,
% drawn in full (every point) by export_cloud_full -- the single renderer used for
% all manuscript clouds. It makes its own square 620 px figure, so the shared `fig`
% used by the 2D figures is left untouched.
mat_1s0_1fs = fullfile(load_base_path, '1s0_1fs', '1s0_1fs.mat');
export_cloud_full(mat_1s0_1fs, '1s0_1fs', 20000, -25, 10, ...
                  fullfile(export_base_path, '1s0_1fs-cloud_partial_200as.pdf'));
export_cloud_full(mat_1s0_1fs, '1s0_1fs', 50000, -25, 10, ...
                  fullfile(export_base_path, '1s0_1fs-cloud_partial_500as.pdf'));
figure(fig);   % restore the shared figure as current for the 2D figures below

%%%%%%%%%%%%%%%%%%%%%%
load(fullfile(load_base_path, '1s0_1ps/1s0_1ps.mat'));
recompute_analytic_distributions;

clf;
i=500;
pdf_name='1s0_1ps-radial_hist_500fs.pdf';
plot_radial_distribution(hist_bins_r, hist_counts_r_traj(i,:), r_values, P_analytic_r);
bump_dist_fonts(20);
exportgraphics(fig, fullfile(export_base_path, pdf_name), 'ContentType', 'vector');
fprintf('%s -> t=%.1f fs\n', pdf_name, i*frame_step*dt/1e-15);

clf;
i=500;
pdf_name='1s0_1ps-polar_hist_500fs.pdf';
plot_polar_distribution(hist_bins_theta, hist_counts_theta_traj(i,:), theta_values, P_theta_analytic);
bump_dist_fonts(20);
exportgraphics(fig, fullfile(export_base_path, pdf_name), 'ContentType', 'vector');
fprintf('%s -> t=%.1f fs\n', pdf_name, i*frame_step*dt/1e-15);

clf;
i=500;
pdf_name='1s0_1ps-azimuthal_hist_500fs.pdf';
plot_azimuthal_distribution(hist_bins_phi, hist_counts_phi_traj(i,:), phi_values, P_phi_analytic);
bump_dist_fonts(20);
exportgraphics(fig, fullfile(export_base_path, pdf_name), 'ContentType', 'vector');
fprintf('%s -> t=%.1f fs\n', pdf_name, i*frame_step*dt/1e-15);

% Manuscript Fig. 2a: the 1s cloud over the first 50 fs, drawn in full (all 5e6
% points, no decimation) to match the full 2p0 cloud beside it (Fig. 2b).
export_cloud_full(fullfile(load_base_path, '1s0_1ps', '1s0_1ps.mat'), '1s0_1ps', ...
                  5000000, -25, 5, fullfile(export_base_path, '1s0_1ps-cloud_50fs.pdf'));
figure(fig);   % restore the shared figure as current for the 2D figures below

%%%%%%%%%%%%%%%%%%%%%%
load(fullfile(load_base_path, '2p0_10ps/2p0_10ps.mat'));
recompute_analytic_distributions;

clf;
i=500;
pdf_name='2p0_10ps-radial_hist_5ps.pdf';
plot_radial_distribution(hist_bins_r, hist_counts_r_traj(i,:), r_values, P_analytic_r);
bump_dist_fonts(20);
exportgraphics(fig, fullfile(export_base_path, pdf_name), 'ContentType', 'vector');
fprintf('%s -> t=%.1f fs\n', pdf_name, i*frame_step*dt/1e-15);

clf;
i=500;
pdf_name='2p0_10ps-polar_hist_5ps.pdf';
plot_polar_distribution(hist_bins_theta, hist_counts_theta_traj(i,:), theta_values, P_theta_analytic);
bump_dist_fonts(20);
exportgraphics(fig, fullfile(export_base_path, pdf_name), 'ContentType', 'vector');
fprintf('%s -> t=%.1f fs\n', pdf_name, i*frame_step*dt/1e-15);

clf;
i=500;
pdf_name='2p0_10ps-azimuthal_hist_5ps.pdf';
plot_azimuthal_distribution(hist_bins_phi, hist_counts_phi_traj(i,:), phi_values, P_phi_analytic);
bump_dist_fonts(20);
exportgraphics(fig, fullfile(export_base_path, pdf_name), 'ContentType', 'vector');
fprintf('%s -> t=%.1f fs\n', pdf_name, i*frame_step*dt/1e-15);

clf;
pdf_name='2p0_10ps-distribution_deviation.pdf';
plot_distributions_deviation(time_arr, dev_arr_R, dev_arr_theta, dev_arr_phi);
bump_dist_fonts(18);   % Fig 8 panels are 0.40\textwidth -> slightly smaller font
exportgraphics(fig, fullfile(export_base_path, pdf_name), 'ContentType', 'vector');
fprintf('%s -> t=%.1f fs\n', pdf_name, time_arr(1,end)/1e-15);

clf;
% Manuscript Fig. 10 (left): 2p0 energies over the first 1 ps of the 2p0_10ps
% run. The filename must match the manuscript \includegraphics
% (2p0_1ps-energies.pdf); do not reuse the 2p_m0_1ps-energies.pdf name, which
% is written further below from the separate 2p_m0_1ps run.
pdf_name='2p0_1ps-energies.pdf';
idx=1:size(time_arr,2)/10;
plot_energies_with_bands(n, l, m, time_arr(idx), ...
    KE_radial_traj_M(:, idx), KE_theta_traj_M(:, idx), KE_phi_traj_M(:, idx), ...
    V_traj_M(:, idx), E_traj_M(:, idx));
bump_dist_fonts(18);   % Fig 10 (left) -- 0.40\textwidth panel
exportgraphics(fig, fullfile(export_base_path, pdf_name), 'ContentType', 'vector');
fprintf('%s -> t=%.1f fs\n', pdf_name, time_arr(1,end/10)/1e-15);

clf;
i=10;
pdf_name='2p0_10ps-polar_hist_100fs.pdf';
plot_polar_distribution(hist_bins_theta, hist_counts_theta_traj(i,:), theta_values, P_theta_analytic);
bump_dist_fonts(20);
exportgraphics(fig, fullfile(export_base_path, pdf_name), 'ContentType', 'vector');
fprintf('%s -> t=%.1f fs\n', pdf_name, i*frame_step*dt/1e-15);

clf;
i=100;
pdf_name='2p0_10ps-polar_hist_1ps.pdf';
plot_polar_distribution(hist_bins_theta, hist_counts_theta_traj(i,:), theta_values, P_theta_analytic);
bump_dist_fonts(20);
exportgraphics(fig, fullfile(export_base_path, pdf_name), 'ContentType', 'vector');
fprintf('%s -> t=%.1f fs\n', pdf_name, i*frame_step*dt/1e-15);

clf;
i=1000;
pdf_name='2p0_10ps-polar_hist_10ps.pdf';
plot_polar_distribution(hist_bins_theta, hist_counts_theta_traj(i,:), theta_values, P_theta_analytic);
bump_dist_fonts(20);
exportgraphics(fig, fullfile(export_base_path, pdf_name), 'ContentType', 'vector');
fprintf('%s -> t=%.1f fs\n', pdf_name, i*frame_step*dt/1e-15);

% Manuscript Fig. 2b: the 2p0 dumbbell cloud, drawing ALL 2e8 trajectory points
% (Inf cap = the whole run) so the Brownian texture is fully resolved. This cloud is
% the only 2p0 figure that touches the raw trajectory, so we free the just-loaded
% double arrays here and let export_cloud_full re-read them in single precision --
% that keeps the full render near ~5 GB instead of OOM-killing on a 16 GB machine.
% ~12 min, dominated by rasterizing 2e8 segments.
clear r_traj theta_traj phi_traj;
fprintf('Rendering Fig. 2b (full 2e8-point 2p0 cloud, ~12 min)...\n');
export_cloud_full(fullfile(load_base_path, '2p0_10ps', '2p0_10ps.mat'), '2p0_10ps', ...
                  inf, -25, 5, fullfile(export_base_path, '2p0_10ps-cloud_2ps.pdf'));
figure(fig);   % export_cloud_full used its own figure; make the shared one current again

%%%%%%%%%%%%%%%%%%%%%

% Manuscript Fig. 5: radial, polar, and azimuthal distributions of the
% (2,1,1) state, accumulated over the full run (same style as Figs. 3-4).
load(fullfile(load_base_path, '2p_m1_1ps/2p_m1_1ps.mat'));
recompute_analytic_distributions;

clf;
pdf_name='2p_m1_1ps-radial_hist.pdf';
plot_radial_distribution(hist_bins_r, hist_counts_r, r_values, P_analytic_r);
bump_dist_fonts(20);
exportgraphics(fig, fullfile(export_base_path, pdf_name), 'ContentType', 'vector');
fprintf('%s -> t=%.1f fs (full run)\n', pdf_name, n_steps*dt/1e-15);

clf;
pdf_name='2p_m1_1ps-polar_hist.pdf';
plot_polar_distribution(hist_bins_theta, hist_counts_theta, theta_values, P_theta_analytic);
bump_dist_fonts(20);
exportgraphics(fig, fullfile(export_base_path, pdf_name), 'ContentType', 'vector');
fprintf('%s -> t=%.1f fs (full run)\n', pdf_name, n_steps*dt/1e-15);

clf;
pdf_name='2p_m1_1ps-azimuthal_hist.pdf';
plot_azimuthal_distribution(hist_bins_phi, hist_counts_phi, phi_values, P_phi_analytic);
bump_dist_fonts(20);
exportgraphics(fig, fullfile(export_base_path, pdf_name), 'ContentType', 'vector');
fprintf('%s -> t=%.1f fs (full run)\n', pdf_name, n_steps*dt/1e-15);

%%%%%%%%%%%%%%%%%%%%%

% --- extended-only diagnostics: (2,1,m) per-state energies and xy clouds.
% None are used in the manuscript .tex or README. ---
if extended

load(fullfile(load_base_path, '2p_m0_1ps/2p_m0_1ps.mat'));
recompute_analytic_distributions;

clf;
pdf_name='2p_m0_1ps-energies.pdf';
plot_energies_with_bands(n, l, m, time_arr, ...
    KE_radial_traj_M, KE_theta_traj_M, KE_phi_traj_M, V_traj_M, E_traj_M);
exportgraphics(fig, fullfile(export_base_path, pdf_name), 'ContentType', 'vector');
fprintf('%s -> t=%.1f fs\n', pdf_name, time_arr(1,end)/1e-15);

clf;
i=max(2, min(size(r_traj,2), 2000000));
pdf_name='2p_m0_1ps-cloud_xy.pdf';
plot_trajectory_3d(n, l, r_traj, theta_traj, phi_traj, a_0, i, '2p_m0_1ps', true);
set(gca, 'Toolbar', []);   % strip the interactive axes toolbar from the export
exportgraphics(fig, fullfile(export_base_path, pdf_name), 'ContentType', 'image');
fprintf('%s -> t=%.1f fs\n', pdf_name, i*dt/1e-15);

%%%%%%%%%%%%%%%%%%%%%

% Re-load (2,1,1): the histograms above are written in the default path, but the
% extended energies/cloud below still need this run's full workspace, and the
% 2p_m0 block above has since overwritten it.
load(fullfile(load_base_path, '2p_m1_1ps/2p_m1_1ps.mat'));
recompute_analytic_distributions;

clf;
pdf_name='2p_m1_1ps-energies.pdf';
plot_energies_with_bands(n, l, m, time_arr, ...
    KE_radial_traj_M, KE_theta_traj_M, KE_phi_traj_M, V_traj_M, E_traj_M);
exportgraphics(fig, fullfile(export_base_path, pdf_name), 'ContentType', 'vector');
fprintf('%s -> t=%.1f fs\n', pdf_name, time_arr(1,end)/1e-15);

clf;
i=max(2, min(size(r_traj,2), 2000000));
pdf_name='2p_m1_1ps-cloud_xy.pdf';
plot_trajectory_3d(n, l, r_traj, theta_traj, phi_traj, a_0, i, '2p_m1_1ps', true);
set(gca, 'Toolbar', []);   % strip the interactive axes toolbar from the export
exportgraphics(fig, fullfile(export_base_path, pdf_name), 'ContentType', 'image');
fprintf('%s -> t=%.1f fs\n', pdf_name, i*dt/1e-15);

%%%%%%%%%%%%%%%%%%%%%

load(fullfile(load_base_path, '2p_mn1_1ps/2p_mn1_1ps.mat'));
recompute_analytic_distributions;

clf;
pdf_name='2p_mn1_1ps-energies.pdf';
plot_energies_with_bands(n, l, m, time_arr, ...
    KE_radial_traj_M, KE_theta_traj_M, KE_phi_traj_M, V_traj_M, E_traj_M);
exportgraphics(fig, fullfile(export_base_path, pdf_name), 'ContentType', 'vector');
fprintf('%s -> t=%.1f fs\n', pdf_name, time_arr(1,end)/1e-15);

clf;
i=max(2, min(size(r_traj,2), 2000000));
pdf_name='2p_mn1_1ps-cloud_xy.pdf';
plot_trajectory_3d(n, l, r_traj, theta_traj, phi_traj, a_0, i, '2p_mn1_1ps', true);
set(gca, 'Toolbar', []);   % strip the interactive axes toolbar from the export
exportgraphics(fig, fullfile(export_base_path, pdf_name), 'ContentType', 'image');
fprintf('%s -> t=%.1f fs\n', pdf_name, i*dt/1e-15);

end  % extended diagnostics

%%%%%%%%%%%%%%%%%%%%%

load(fullfile(load_base_path, '2s0_10ps/2s0_10ps.mat'));
recompute_analytic_distributions;

clf;
i=30;
pdf_name='2s0_10ps-radial_hist_300fs.pdf';
plot_radial_distribution(hist_bins_r, hist_counts_r_traj(i,:), r_values, P_analytic_r);
bump_dist_fonts(20);
exportgraphics(fig, fullfile(export_base_path, pdf_name), 'ContentType', 'vector');
fprintf('%s -> t=%.1f fs\n', pdf_name, i*frame_step*dt/1e-15);


clf;
i=35;
pdf_name='2s0_10ps-radial_hist_350fs.pdf';
plot_radial_distribution(hist_bins_r, hist_counts_r_traj(i,:), r_values, P_analytic_r);
bump_dist_fonts(20);
exportgraphics(fig, fullfile(export_base_path, pdf_name), 'ContentType', 'vector');
fprintf('%s -> t=%.1f fs\n', pdf_name, i*frame_step*dt/1e-15);

clf;
i=500;
pdf_name='2s0_10ps-radial_hist_5ps.pdf';
plot_radial_distribution(hist_bins_r, hist_counts_r_traj(i,:), r_values, P_analytic_r);
bump_dist_fonts(20);
exportgraphics(fig, fullfile(export_base_path, pdf_name), 'ContentType', 'vector');
fprintf('%s -> t=%.1f fs\n', pdf_name, i*frame_step*dt/1e-15);

clf;
pdf_name='2s0_10ps-distribution_deviation.pdf';
plot_distributions_deviation(time_arr, dev_arr_R, dev_arr_theta, dev_arr_phi);
bump_dist_fonts(18);   % Fig 8 panels are 0.40\textwidth -> slightly smaller font
exportgraphics(fig, fullfile(export_base_path, pdf_name), 'ContentType', 'vector');
fprintf('%s -> t=%.1f fs\n', pdf_name, time_arr(1,end)/1e-15);

clf;
pdf_name='2s0_10ps-energies.pdf';
plot_energies_with_bands(n, l, m, time_arr, ...
    KE_radial_traj_M, KE_theta_traj_M, KE_phi_traj_M, V_traj_M, E_traj_M);
bump_dist_fonts(18);   % Fig 10 (right) -- 0.40\textwidth panel
exportgraphics(fig, fullfile(export_base_path, pdf_name), 'ContentType', 'vector');
fprintf('%s -> t=%.1f fs\n', pdf_name, time_arr(1,end)/1e-15);

function bump_dist_fonts(fs)
    % Enlarge the axis/label/legend fonts of a 2D distribution figure so the small
    % 0.32\textwidth manuscript panels (Figs 3-7) read at body-text size. Applied
    % here, after the plot call -- the shared plot_*_distribution functions (also
    % used by the dashboard movie) keep their default smaller fonts.
    if nargin < 1 || isempty(fs)
        fs = 20;
    end
    set(gca, 'FontSize', fs);
    set([get(gca, 'XLabel'), get(gca, 'YLabel')], 'FontSize', fs + 2);
    lg = get(gca, 'Legend');
    if ~isempty(lg)
        lg.FontSize = fs - 2;
    end
end

function recompute_analytic_distributions
    % Local functions have their own workspace, so pull the per-run state set
    % by the most recent load(...) from the caller (the script) before use.
    r_values     = evalin('caller', 'r_values');
    theta_values = evalin('caller', 'theta_values');
    phi_values   = evalin('caller', 'phi_values');
    n = evalin('caller', 'n');
    l = evalin('caller', 'l');
    m = evalin('caller', 'm');

    R_analytic = radial_wavefunction(r_values, n, l);
    P_analytic_r = r_values.^2 .* abs(R_analytic).^2;

    Y_theta_analytic = angular_wavefunction_theta(theta_values, l, m);
    P_theta_analytic = 2 * pi * sin(theta_values) .* abs(Y_theta_analytic).^2;

    P_phi_analytic = ones(size(phi_values)) / (2 * pi);

    assignin('caller', 'R_analytic', R_analytic);
    assignin('caller', 'P_analytic_r', P_analytic_r);
    assignin('caller', 'Y_theta_analytic', Y_theta_analytic);
    assignin('caller', 'P_theta_analytic', P_theta_analytic);
    assignin('caller', 'P_phi_analytic', P_phi_analytic);
end

function tf = is_scan_parameter_set(params_set_name)
    tf = ~isempty(regexp(params_set_name, '^(2p0|2s0)_scan_vmax_[0-9]+p[0-9]+c$', 'once'));
end

function maybe_analyse_cutoff_scan(params_set_name, export_dir)
    tokens = regexp(params_set_name, '^(2p0|2s0)_scan_vmax_', 'tokens', 'once');
    if isempty(tokens)
        return;
    end
    state_label = tokens{1};
    suffixes = {'0p01c', '0p03c', '0p1c', '0p3c', '1p0c', '2p0c', '3p0c'};

    missing = {};
    for k = 1:length(suffixes)
        run_name = [state_label '_scan_vmax_' suffixes{k}];
        mat_path = fullfile('data', run_name, [run_name '.mat']);
        if ~exist(mat_path, 'file')
            missing{end + 1} = run_name; %#ok<AGROW>
        end
    end

    if isempty(missing)
        analyse_cutoff_scan(state_label, export_dir);
    else
        fprintf('Cutoff scan for %s is not complete yet (%d/%d files present).\n', ...
                state_label, length(suffixes) - length(missing), length(suffixes));
        fprintf('Missing run(s):\n');
        for k = 1:length(missing)
            fprintf('  %s\n', missing{k});
        end
    end
end
