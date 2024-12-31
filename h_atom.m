clearvars -except params_set_name;

rng(7); % to reproduce exactly manuscript figures

if ~exist('params_set_name', 'var')
    params_set_name =  'test'; % default value can be changed
end

[n, l, m, n_steps, dt, traj_points, sigma_r_factor] = parameters(params_set_name);

M = 1;           % number of particles
n_frames = 1000; % Total number of video frames
max_traj=min(n_steps, traj_points); % limit for minimal memory

% Discretization Constants
dr = 1e-15;
dtheta = 1e-5;
dphi = 1e-5;

base_dir = strcat('data/', params_set_name);
status = mkdir(base_dir);
if status
    disp(['Directory "', base_dir, '" created successfully.']);
else
    disp(['Failed to create directory "', base_dir, '".']);
    exit;
end

video_filename = fullfile(base_dir, strcat(params_set_name,'.mp4'));
frame_step = (n_steps / n_frames);

% Energy Variables
S_KE_radial = 0;
S_KE_theta = 0;
S_KE_phi = 0;
S_V = 0;
S_E = 0;

% Physical Constants
[hbar, m_e, a_0, e_charge, epsilon_0, c] = constants();
sigma = sqrt(hbar / m_e);

[r, theta, phi] = initial_coordinates(n, l, m, M, a_0, params_set_name);

time_arr = linspace(0, n_steps * dt, n_frames);

r_traj = zeros(1, max_traj);
theta_traj = zeros(1, max_traj);
phi_traj = zeros(1, max_traj);

KE_radial_traj = zeros(1, n_frames);
KE_theta_traj = zeros(1, n_frames);
KE_phi_traj = zeros(1, n_frames);
V_traj = zeros(1, n_frames);
E_traj = zeros(1, n_frames);

dev_arr_R = zeros(1, n_frames);
dev_arr_theta = zeros(1, n_frames);
dev_arr_phi = zeros(1, n_frames);

% Small values to prevent division by zero
epsilon_r = 1e-4 * a_0;
epsilon_theta = 1e-7;
epsilon_R = 1e-7;
epsilon_Y_theta = 1e-7;
epsilon_Y_phi = 1e-7;

% Histogram Parameters for r
[r_min, r_max] = radial_histogram_range(n, l, a_0, params_set_name);
n_bins_r = 100; 
hist_bins_r = linspace(r_min, r_max, n_bins_r);
hist_counts_r = zeros(1, n_bins_r - 1);
hist_counts_r_traj = zeros(n_frames, n_bins_r - 1);

% Histogram Parameters for theta
theta_min = 0;
theta_max = pi;
n_bins_theta = 100;
hist_bins_theta = linspace(theta_min, theta_max, n_bins_theta);
hist_counts_theta = zeros(1, n_bins_theta - 1);
hist_counts_theta_traj = zeros(n_frames, n_bins_theta - 1);

% Histogram Parameters for phi
phi_min = 0;
phi_max = 2 * pi;
n_bins_phi = 100;
hist_bins_phi = linspace(phi_min, phi_max, n_bins_phi);
hist_counts_phi = zeros(1, n_bins_phi - 1);
hist_counts_phi_traj = zeros(n_frames, n_bins_phi - 1);

% Precompute Analytic Radial Distribution
r_values = linspace(r_min, r_max, 1000);
R_analytic = radial_wavefunction(r_values, n, l);
P_analytic_r = r_values.^2 .* abs(R_analytic).^2;

% Precompute Analytic Angular Distributions
theta_values = linspace(theta_min, theta_max, 1000);
Y_theta_analytic = angular_wavefunction_theta(theta_values, l, m);
P_theta_analytic = 2 * pi * sin(theta_values) .* abs(Y_theta_analytic).^2;

% Precompute Analytic Azimuthal Distribution
phi_values = linspace(phi_min, phi_max, 1000);
Y_phi_analytic = angular_wavefunction_phi(phi_values, m);
P_phi_analytic = abs(Y_phi_analytic).^2 / (2 * pi);

% Initialize Video Writer
v = VideoWriter(video_filename, 'MPEG-4');
v.FrameRate = 10;
open(v);

frame_index = 1;
startTime = datetime('now');
fprintf('Start time: %s\n', startTime);

fig = figure;
        
tic;
for i = 1:n_steps
    % Generate Gaussian Increments for Radial Component (r) of Brownian Motion
    dW_r = randn(M, 1) * sigma_r_factor *sigma * sqrt(dt);

    % Compute Radial Wave Function and Its Derivative at r
    R = radial_wavefunction(r, n, l);
    dR_dr = derivative(@(r) radial_wavefunction(r, n, l), r, dr);

    % Prevent Division by Zero or Very Small r Values
    r_safe = max(r, epsilon_r);

    % Handle R_safe to Prevent Division by Zero
    R_safe = max(abs(R), epsilon_R) .* (1 - 2 * (R < 0)); % sign(R);

    % Compute Radial Drift Velocity
    v_r_drift = (hbar / m_e) * (dR_dr ./ R_safe);
    v_r_curv = (hbar / m_e) * (1 ./ r_safe);
    
    % Limit the drift velocity to ensure polar energy convergence
    v_r_drift_idx = find(abs(v_r_drift) > c);
    v_r_drift(v_r_drift_idx) = c * sign(v_r_drift(v_r_drift_idx));

    % Limit the drift velocity (due to curviture) to ensure radial energy convergence
    v_r_curv_idx = find(abs(v_r_curv) > c);
    v_r_curv(v_r_curv_idx) = c * sign(v_r_curv(v_r_curv_idx));

    % Update the Radial Position
    r = r + (v_r_curv + v_r_drift) .* dt + dW_r;

    % Reflective Boundary at r = 0
    r = abs(r);

    % Accumulate Histogram Counts for R
    hist_counts_r = hist_counts_r + histcounts(r, hist_bins_r);

    % Generate Gaussian Increments for Brownian Motion
    dW_theta = (sigma * sqrt(dt)) .* randn(M, 1);

    % Prevent theta from Reaching Exactly 0 or pi
    theta_safe = max(min(theta, pi - epsilon_theta), epsilon_theta);

    % Compute Polar Wavefunction and Its Derivative at theta
    Y_theta = angular_wavefunction_theta(theta_safe, l, m);
    dY_theta_dtheta = derivative(@(theta) angular_wavefunction_theta(theta, l, m), theta_safe, dtheta);

    % Handle Y_theta_safe to Prevent Division by Zero
    Y_theta_safe = max(abs(Y_theta), epsilon_Y_theta) .* (1 - 2 * (Y_theta < 0)); % sign(Y_theta);

    % Compute Polar Drift Velocity
    v_theta_drift = (hbar / m_e) * (dY_theta_dtheta ./ Y_theta_safe) ./ r_safe;

    % Limit the drift velocity to ensure polar energy convergence
    v_theta_drift_idx = find(abs(v_theta_drift) > c);
    v_theta_drift(v_theta_drift_idx) = c * sign(v_theta_drift(v_theta_drift_idx));

    v_theta_curv = (hbar / m_e) * (cot(theta_safe) / 2 ) ./ r_safe;

    % Limit the drift velocity (due to curviture) to ensure polar energy convergence
    v_theta_curv_idx = find(abs(v_theta_curv) > c);
    v_theta_curv(v_theta_curv_idx) = c * sign(v_theta_curv(v_theta_curv_idx));
   
    % Update the Polar Angle
    theta = theta + (1./ r_safe) .* ((v_theta_curv + v_theta_drift) .* dt + dW_theta);
  
    % Accumulate Histogram Counts for theta
    hist_counts_theta = hist_counts_theta + histcounts(theta,  hist_bins_theta);
    
    % Generate Random Steps (Gaussian Increments for Brownian Motion)
    dW_phi = (sigma * sqrt(dt)) .* randn(M, 1);

    % Compute Azimuthal Wavefunction and Its Derivative at phi
    Y_phi = angular_wavefunction_phi(phi, m);
    dY_phi_dphi = derivative(@(phi) angular_wavefunction_phi(phi, m), phi, dphi);

    % Handle Y_phi_safe to Prevent Division by Zero
    Y_phi_safe = max(abs(Y_phi), epsilon_Y_phi) .* (1 - 2 * (Y_phi < 0)); % sign(Y_phi);

    % Compute Azimuthal Drift Velocity
    v_phi_drift = (hbar / m_e) * (dY_phi_dphi ./ Y_phi_safe) ./ (r_safe .* sin(theta_safe));

    % Limit the drift velocity to ensure polar energy convergence
    v_phi_drift_idx = find(abs(v_phi_drift) > c);
    v_phi_drift(v_phi_drift_idx) = c * sign(v_phi_drift(v_phi_drift_idx));

    % Update Azimuthal Angle
    phi = phi + 1./(r_safe .* sin(theta_safe)) .* (v_phi_drift .* dt + dW_phi);

    % Adjust phi to be Within [0, 2*pi]
    phi = mod(phi, 2 * pi);

    % Accumulate Histogram Counts for phi
    hist_counts_phi = hist_counts_phi + histcounts(phi, hist_bins_phi);

    % store particle coordinates
    if i <= max_traj && M == 1
        r_traj(i) = r;
        theta_traj(i) = theta;
        phi_traj(i) = phi;
    end

    % Compute Kinetic Energy Using Drift Velocities
    KE_radial = 0.5 * m_e *  (v_r_drift.^2); 
    KE_theta = 0.5 * m_e * (v_theta_drift.^2); 
    KE_phi = 0.5 * m_e * (v_phi_drift.^2); 

    % Compute Potential Energy
    V = -e_charge^2 ./ (4 * pi * epsilon_0 * r_safe);
    
    % Compute Total Energy
    E = KE_radial + KE_theta + KE_phi + V;
    
    % Accumulate Energies
    S_KE_radial = S_KE_radial + sum(KE_radial);
    S_KE_theta = S_KE_theta + sum(KE_theta);
    S_KE_phi = S_KE_phi + sum(KE_phi);
    S_V = S_V + sum(V);
    S_E = S_E + sum(E);

    if mod(i, frame_step) == 0
        fprintf('i=%d, t=%.5e, KE_r=%.5e, K_azim=%.5e, V=%.5e, E=%.5e\r', ...
            i, i*dt, S_KE_radial/i, S_KE_theta / i, S_V/i, S_E/i);

        KE_radial_traj(frame_index) = S_KE_radial /  i;
        KE_theta_traj(frame_index) = S_KE_theta /  i;
        KE_phi_traj(frame_index) = S_KE_phi / i;
        V_traj(frame_index) = S_V / i;
        E_traj(frame_index) = S_E / i;
        hist_counts_r_traj(frame_index, :) = hist_counts_r; 
        hist_counts_theta_traj(frame_index, :) = hist_counts_theta;
        hist_counts_phi_traj(frame_index, :) = hist_counts_phi;
    
        dev_arr_R(frame_index) = compute_distributions_deviation(...
            hist_bins_r, hist_counts_r, r_values, P_analytic_r);

        dev_arr_theta(frame_index) = compute_distributions_deviation(...
            hist_bins_theta, hist_counts_theta, theta_values, P_theta_analytic);

        dev_arr_phi(frame_index) = compute_distributions_deviation(...
                hist_bins_phi, hist_counts_phi, phi_values, P_phi_analytic);        
        
        clf;
        set(gcf, 'Position', [50, 50, 1600, 900]);
        subplot(2,3,1);
        plot_radial_distribution(hist_bins_r, hist_counts_r, r_values, P_analytic_r);
        subplot(2,3,2);
        plot_polar_distribution(hist_bins_theta, hist_counts_theta, theta_values, P_theta_analytic);
        subplot(2,3,3);
        plot_azimuthal_distribution(hist_bins_phi, hist_counts_phi, phi_values, P_phi_analytic);
        subplot(2,3,4);

        idx = 1:frame_index;
        plot_distributions_deviation(time_arr(idx), dev_arr_R(idx), dev_arr_theta(idx), dev_arr_phi(idx));
        subplot(2,3,5);
        plot_energies(n, l, m, time_arr(idx), KE_radial_traj(idx), ...
            KE_theta_traj(idx), KE_phi_traj(idx), V_traj(idx), E_traj(idx));

        subplot(2,3,6);
        plot_trajectory_3d(n, l, r_traj, theta_traj, phi_traj, a_0, i, params_set_name);

        annotation_text = sprintf('t = %.2e s', time_arr(frame_index));
        annotation('textbox', [0.02, 0.95, 0.1, 0.05], 'String', annotation_text, ...
        'Units', 'normalized', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', ...
        'FontSize', 12, 'FontWeight', 'bold', 'EdgeColor', 'none');

        frame = getframe(gcf);
        writeVideo(v, frame);
        frame_index = frame_index + 1;
    end
end

totalElapsedTime = toc;
close(v);
fprintf('Total elapsed time: %.4f seconds\n', totalElapsedTime);

exportgraphics(fig, fullfile(base_dir, strcat(params_set_name,'.pdf')), 'ContentType', 'image');
clear fig;

save(fullfile(base_dir, strcat(params_set_name,'.mat')), '-v7.3');
