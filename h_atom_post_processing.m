clear;

export_base_path = 'figures';
load_base_path = 'data';

status = mkdir(export_base_path);

if status
    disp(['Directory "', export_base_path, '" created successfully.']);
else
    disp(['Failed to create directory "', export_base_path, '".']);
    exit;
end

fig = figure;
defaultFontSize = 20;

% Set default font properties for the figure
set(fig, 'DefaultAxesFontSize', defaultFontSize, ...
         'DefaultTextFontSize', defaultFontSize, ...
         'DefaultLegendFontSize', defaultFontSize);

%%%%%%%%%%%%%%%%%%%%%

load(fullfile(load_base_path, '1s0_1fs/1s0_1fs.mat')); % fig_1

clf;
i=20000;
pdf_name='1s0_1fs-cloud_partial_200as.pdf';
plot_trajectory_3d(n, l, r_traj, theta_traj, phi_traj, a_0, i, params_set_name, true);
view(-25, 10);
exportgraphics(fig, fullfile(export_base_path, pdf_name), 'ContentType', 'image');
fprintf('%s -> t=%.1f fs\n', pdf_name, i*dt/1e-15);

clf;
i=50000;
pdf_name='1s0_1fs-cloud_partial_500as.pdf';
plot_trajectory_3d(n, l, r_traj, theta_traj, phi_traj, a_0, i, params_set_name, true);
view(-25, 10);
exportgraphics(fig, fullfile(export_base_path, pdf_name), 'ContentType', 'image');
fprintf('%s -> t=%.1f fs\n', pdf_name, i*dt/1e-15);

%%%%%%%%%%%%%%%%%%%%%%
load(fullfile(load_base_path, '1s0_1ps/1s0_1ps.mat'));

clf;
i=500;
pdf_name='1s0_1ps-radial_hist_500fs.pdf';
plot_radial_distribution(hist_bins_r, hist_counts_r_traj(i,:), r_values, P_analytic_r);
exportgraphics(fig, fullfile(export_base_path, pdf_name), 'ContentType', 'vector');
fprintf('%s -> t=%.1f fs\n', pdf_name, i*frame_step*dt/1e-15);

clf;
i=500;
pdf_name='1s0_1ps-polar_hist_500fs.pdf';
plot_polar_distribution(hist_bins_theta, hist_counts_theta_traj(i,:), theta_values, P_theta_analytic);
exportgraphics(fig, fullfile(export_base_path, pdf_name), 'ContentType', 'vector');
fprintf('%s -> t=%.1f fs\n', pdf_name, i*frame_step*dt/1e-15);

clf;
i=500;
pdf_name='1s0_1ps-azimuthal_hist_500fs.pdf';
plot_azimuthal_distribution(hist_bins_phi, hist_counts_phi_traj(i,:), phi_values, P_phi_analytic);
exportgraphics(fig, fullfile(export_base_path, pdf_name), 'ContentType', 'vector');
fprintf('%s -> t=%.1f fs\n', pdf_name, i*frame_step*dt/1e-15);

clf;
i=5000000;
pdf_name='1s0_1ps-cloud_50fs.pdf';
plot_trajectory_3d(n, l, r_traj, theta_traj, phi_traj, a_0, i, params_set_name, true);
exportgraphics(fig, fullfile(export_base_path, pdf_name), 'ContentType', 'image');
fprintf('%s -> t=%.1f fs\n', pdf_name, i*dt/1e-15);

%%%%%%%%%%%%%%%%%%%%%%
load(fullfile(load_base_path, '2p0_10ps/2p0_10ps.mat'));

clf;
i=500;
pdf_name='2p0_10ps-radial_hist_5ps.pdf';
plot_radial_distribution(hist_bins_r, hist_counts_r_traj(i,:), r_values, P_analytic_r);
exportgraphics(fig, fullfile(export_base_path, pdf_name), 'ContentType', 'vector');
fprintf('%s -> t=%.1f fs\n', pdf_name, i*frame_step*dt/1e-15);

clf;
i=500;
pdf_name='2p0_10ps-polar_hist_5ps.pdf';
plot_polar_distribution(hist_bins_theta, hist_counts_theta_traj(i,:), theta_values, P_theta_analytic);
exportgraphics(fig, fullfile(export_base_path, pdf_name), 'ContentType', 'vector');
fprintf('%s -> t=%.1f fs\n', pdf_name, i*frame_step*dt/1e-15);

clf;
i=500;
pdf_name='2p0_10ps-azimuthal_hist_5ps.pdf';
plot_azimuthal_distribution(hist_bins_phi, hist_counts_phi_traj(i,:), phi_values, P_phi_analytic);
exportgraphics(fig, fullfile(export_base_path, pdf_name), 'ContentType', 'vector');
fprintf('%s -> t=%.1f fs\n', pdf_name, i*frame_step*dt/1e-15);

clf;
pdf_name='2p0_10ps-distribution_deviation.pdf';
plot_distributions_deviation(time_arr, dev_arr_R, dev_arr_theta, dev_arr_phi);
exportgraphics(fig, fullfile(export_base_path, pdf_name), 'ContentType', 'vector');
fprintf('%s -> t=%.1f fs\n', pdf_name, time_arr(1,end)/1e-15);

clf;
pdf_name='2p0_1ps-energies.pdf';
idx=1:size(time_arr,2)/10;
plot_energies(n, l, m, time_arr(idx), KE_radial_traj(idx), ...
    KE_theta_traj(idx), KE_phi_traj(idx), V_traj(idx), E_traj(idx));
exportgraphics(fig, fullfile(export_base_path, pdf_name), 'ContentType', 'vector');
fprintf('%s -> t=%.1f fs\n', pdf_name, time_arr(1,end/10)/1e-15);

clf;
i=10;
pdf_name='2p0_10ps-polar_hist_100fs.pdf';
plot_polar_distribution(hist_bins_theta, hist_counts_theta_traj(i,:), theta_values, P_theta_analytic);
exportgraphics(fig, fullfile(export_base_path, pdf_name), 'ContentType', 'vector');
fprintf('%s -> t=%.1f fs\n', pdf_name, i*frame_step*dt/1e-15);

clf;
i=100;
pdf_name='2p0_10ps-polar_hist_1ps.pdf';
plot_polar_distribution(hist_bins_theta, hist_counts_theta_traj(i,:), theta_values, P_theta_analytic);
exportgraphics(fig, fullfile(export_base_path, pdf_name), 'ContentType', 'vector');
fprintf('%s -> t=%.1f fs\n', pdf_name, i*frame_step*dt/1e-15);

clf;
i=1000;
pdf_name='2p0_10ps-polar_hist_10ps.pdf';
plot_polar_distribution(hist_bins_theta, hist_counts_theta_traj(i,:), theta_values, P_theta_analytic);
exportgraphics(fig, fullfile(export_base_path, pdf_name), 'ContentType', 'vector');
fprintf('%s -> t=%.1f fs\n', pdf_name, i*frame_step*dt/1e-15);

clf;
i=200000000;
fprintf('Ploting 2p state might took 2-3 min....\n')
pdf_name='2p0_10ps-cloud_2ps.pdf';
plot_trajectory_3d(n, l, r_traj, theta_traj, phi_traj, a_0, i, '2p0_10ps', true);
view(-25, 5);
fprintf('Exporting 2p state might took 2-3 min....\n')
exportgraphics(fig, fullfile(export_base_path, pdf_name), 'ContentType', 'image');
fprintf('%s -> t=%.1f fs\n', pdf_name, i*dt/1e-15);

%%%%%%%%%%%%%%%%%%%%%

load(fullfile(load_base_path, '2s0_10ps/2s0_10ps.mat'));

clf;
i=30;
pdf_name='2s0_10ps-radial_hist_300fs.pdf';
plot_radial_distribution(hist_bins_r, hist_counts_r_traj(i,:), r_values, P_analytic_r);
exportgraphics(fig, fullfile(export_base_path, pdf_name), 'ContentType', 'vector');
fprintf('%s -> t=%.1f fs\n', pdf_name, i*frame_step*dt/1e-15);


clf;
i=35;
pdf_name='2s0_10ps-radial_hist_350fs.pdf';
plot_radial_distribution(hist_bins_r, hist_counts_r_traj(i,:), r_values, P_analytic_r);
exportgraphics(fig, fullfile(export_base_path, pdf_name), 'ContentType', 'vector');
fprintf('%s -> t=%.1f fs\n', pdf_name, i*frame_step*dt/1e-15);

clf;
i=500;
pdf_name='2s0_10ps-radial_hist_5ps.pdf';
plot_radial_distribution(hist_bins_r, hist_counts_r_traj(i,:), r_values, P_analytic_r);
exportgraphics(fig, fullfile(export_base_path, pdf_name), 'ContentType', 'vector');
fprintf('%s -> t=%.1f fs\n', pdf_name, i*frame_step*dt/1e-15);

clf;
pdf_name='2s0_10ps-distribution_deviation.pdf';
plot_distributions_deviation(time_arr, dev_arr_R, dev_arr_theta, dev_arr_phi);
exportgraphics(fig, fullfile(export_base_path, pdf_name), 'ContentType', 'vector');
fprintf('%s -> t=%.1f fs\n', pdf_name, time_arr(1,end)/1e-15);

clf;
pdf_name='2s0_10ps-energies.pdf';
plot_energies(n, l, m, time_arr, KE_radial_traj, KE_theta_traj, KE_phi_traj, V_traj, E_traj);
exportgraphics(fig, fullfile(export_base_path, pdf_name), 'ContentType', 'vector');
fprintf('%s -> t=%.1f fs\n', pdf_name, time_arr(1,end)/1e-15);
