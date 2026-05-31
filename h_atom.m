resume_requested = exist('resume_from_workspace', 'var') && resume_from_workspace;
if resume_requested
    can_resume_same_run = exist('active_params_set_name', 'var') && ...
        exist('params_set_name', 'var') && strcmp(params_set_name, active_params_set_name);
    already_finished = exist('i', 'var') && exist('n_steps', 'var') && i >= n_steps;
    resume_requested = can_resume_same_run && ~already_finished;
end
resume_from_workspace = false;
if ~resume_requested
clearvars -except params_set_name resume_requested resume_from_workspace independent_random_design random_seed dt_override;
resume_from_workspace = false;

if ~exist('independent_random_design', 'var')
    independent_random_design = false;
end
if ~exist('random_seed', 'var') || isempty(random_seed)
    random_seed = 7;
end
if independent_random_design
    rng(random_seed);
    random_design_label = sprintf('independent_seed_%d', random_seed);
else
    random_seed = 7;
    rng(random_seed); % common-random-number design for reproducible sweeps
    random_design_label = 'common_seed_7';
end

if ~exist('params_set_name', 'var')
    params_set_name =  'test'; % default value can be changed
end
active_params_set_name = params_set_name;

[n, l, m, n_steps, dt, traj_points, sigma_r_factor, M, v_max_over_c] = parameters(params_set_name);

% Optional dt sweep: set dt_override in the workspace before running to use a
% different time step than the one defined in parameters.m. n_steps is rescaled
% so the total simulated time (n_steps * dt) is held fixed -- this keeps dt
% sweeps at equal physical time. Clear dt_override to fall back to parameters.m.
if exist('dt_override', 'var') && ~isempty(dt_override)
    total_time = n_steps * dt;
    dt = dt_override;
    n_steps = round(total_time / dt);
    traj_points = min(traj_points, n_steps);
    fprintf(['dt override active: dt = %.3e s, n_steps rescaled to %d ' ...
             '(total time %.3e s held fixed)\n'], dt, n_steps, total_time);
end

% M is supplied by parameters.m so scan runs can use fewer trajectories
% without changing the production presets.
is_cutoff_scan_run = ~isempty(regexp(params_set_name, '_scan_vmax_', 'once'));
is_dt_scan_run = ~isempty(regexp(params_set_name, '_scan_dt_', 'once'));
% Both sweeps use the same kinetic-error-vs-reference plotting layout.
is_kinetic_scan_run = is_cutoff_scan_run || is_dt_scan_run;
output_run_name = params_set_name;
if independent_random_design && is_cutoff_scan_run
    scan_tokens = regexp(params_set_name, '^(2p0|2s0)_scan_vmax_([0-9]+p[0-9]+c)$', ...
                         'tokens', 'once');
    if isempty(scan_tokens)
        error('Independent random design is only defined for cutoff scan runs.');
    end
    output_run_name = sprintf('%s_scan_independent_seed_%d_vmax_%s', ...
                              scan_tokens{1}, random_seed, scan_tokens{2});
elseif independent_random_design
    output_run_name = sprintf('%s_independent_seed_%d', params_set_name, random_seed);
end
% Tag the output folder/files with the dt value so each dt sweep point is saved
% separately instead of overwriting the baseline run's data directory. The dt
% is expressed in zeptoseconds (1 zs = 1e-21 s) for a readable tag, with 'p' in
% place of any decimal point, e.g. dt=1e-20 s -> '_dt_10zs', 5e-22 s -> '_dt_0p5zs'.
if exist('dt_override', 'var') && ~isempty(dt_override)
    dt_zs_str = strrep(sprintf('%g', dt / 1e-21), '.', 'p');
    output_run_name = sprintf('%s_dt_%szs', output_run_name, dt_zs_str);
end
make_live_plots = true;
% By default the first 8 subplots (radial/polar/azimuthal histograms,
% trajectory views, distribution deviation, energies-vs-time, azimuthal
% current) use a single trajectory -- matching the original M = 1
% behavior. The row-3 per-trajectory energy histograms and the end-of-run
% diagnostic always use all M trajectories. Set to false to pool all M
% particles into the main plots as well (cross-particle means + SEM bands).
single_traj_main_plots = true;
n_frames = 1000; % Total number of video frames
if strcmp(params_set_name, '2p_m1_1ps') || strcmp(params_set_name, '2p1_1ps') || strcmp(params_set_name, '2p_m0_1ps') || strcmp(params_set_name, '2p0_1ps') || strcmp(params_set_name, '2p_mn1_1ps')
    n_frames = 80;
end
max_traj=min(n_steps, traj_points); % limit for minimal memory

% Discretization Constants
dr = 1e-15;
dtheta = 1e-5;

base_dir = fullfile('data', output_run_name);
if ~exist(base_dir, 'dir')
    status = mkdir(base_dir);
    if status
        disp(['Directory "', base_dir, '" created successfully.']);
    else
        disp(['Failed to create directory "', base_dir, '".']);
        exit;
    end
else
    disp(['Directory "', base_dir, '" already exists. Existing outputs may be overwritten.']);
end

video_filename = fullfile(base_dir, strcat(output_run_name,'.mp4'));
frame_step = (n_steps / n_frames);
fprintf('Random design: %s (seed = %d)\n', random_design_label, random_seed);

% Energy Variables
S_KE_radial = 0;
S_KE_theta = 0;
S_KE_phi = 0;
S_b_phi = 0;
S_V = 0;
S_E = 0;

% Per-particle energy accumulators (size [M,1]) for SEM across trajectories
S_KE_radial_M = zeros(M, 1);
S_KE_theta_M  = zeros(M, 1);
S_KE_phi_M    = zeros(M, 1);
S_V_M         = zeros(M, 1);
S_E_M         = zeros(M, 1);

% Physical Constants
[hbar, m_e, a_0, e_charge, epsilon_0, c] = constants();
sigma = sqrt(hbar / m_e);
v_max = v_max_over_c * c;

[r, theta, phi] = initial_coordinates(n, l, m, M, a_0, params_set_name);
phi_unwrapped = phi;

% Cutoff-scan diagnostics. The nodal-set diagnostics are enabled only for
% the two states used in Appendix E; cutoff engagement is recorded for all
% states.
N_cross_per_traj = zeros(M, 1);
occupation_count_per_traj = zeros(M, 1);
cutoff_engagement_r_count_per_traj = zeros(M, 1);
cutoff_engagement_theta_count_per_traj = zeros(M, 1);
cutoff_engagement_phi_count_per_traj = zeros(M, 1);

use_polar_node_diagnostic = (n == 2 && l == 1 && m == 0);
use_radial_node_diagnostic = (n == 2 && l == 0 && m == 0);
if use_polar_node_diagnostic
    previous_node_value = cos(theta);
elseif use_radial_node_diagnostic
    previous_node_value = r - 2 * a_0;
else
    previous_node_value = NaN(M, 1);
end

time_arr = (1:n_frames) * frame_step * dt;

r_traj = zeros(1, max_traj);
theta_traj = zeros(1, max_traj);
phi_traj = zeros(1, max_traj);

KE_radial_traj = zeros(1, n_frames);
KE_theta_traj = zeros(1, n_frames);
KE_phi_traj = zeros(1, n_frames);
b_phi_avg_traj = zeros(1, n_frames);
V_traj = zeros(1, n_frames);
E_traj = zeros(1, n_frames);

% Per-trajectory running-average time series (size [M, n_frames]).
% Used by plot_energies_with_bands.m to draw +/- SEM bands.
KE_radial_traj_M = zeros(M, n_frames);
KE_theta_traj_M  = zeros(M, n_frames);
KE_phi_traj_M    = zeros(M, n_frames);
V_traj_M         = zeros(M, n_frames);
E_traj_M         = zeros(M, n_frames);

dev_arr_R = zeros(1, n_frames);
dev_arr_theta = zeros(1, n_frames);
dev_arr_phi = zeros(1, n_frames);

% Small values to prevent division by zero
epsilon_r = 1e-4 * a_0;
epsilon_theta = 1e-7;
epsilon_R = 1e-7;
epsilon_Y_theta = 1e-7;

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
P_phi_analytic = ones(size(phi_values)) / (2 * pi);

if make_live_plots
    % Initialize Video Writer
    v = VideoWriter(video_filename, 'MPEG-4');
    v.FrameRate = 10;
    if strcmp(params_set_name, '2p_m1_1ps') || strcmp(params_set_name, '2p1_1ps') || strcmp(params_set_name, '2p_m0_1ps') || strcmp(params_set_name, '2p0_1ps') || strcmp(params_set_name, '2p_mn1_1ps')
        v.FrameRate = 24;
    end
    open(v);
    fig = figure('Units', 'pixels', 'Position', [50, 50, 2100, 1350]);
else
    v = [];
    fig = [];
end

frame_index = 1;
startTime = datetime('now');
fprintf('Start time: %s\n', startTime);

% Analytical references for the per-trajectory energy histograms.
ref_energies = energy_reference_values(n, l, m);
if n == 2 && l == 1 && m == 0
    kinetic_scan_label = '$\left< T_\theta \right>$';
    kinetic_scan_ref = ref_energies.E_theta;
elseif n == 2 && l == 0 && m == 0
    kinetic_scan_label = '$\left< T_r \right>$';
    kinetic_scan_ref = ref_energies.E_r;
else
    kinetic_scan_label = '$\left< E \right>$';
    kinetic_scan_ref = ref_energies.E;
end

kinetic_mean_relerr_traj = zeros(1, n_frames);
kinetic_median_relerr_traj = zeros(1, n_frames);
kinetic_sem_norm_traj = zeros(1, n_frames);
kinetic_iqr_norm_traj = zeros(1, n_frames);
N_cross_total_traj = zeros(1, n_frames);
N_cross_mean_traj = zeros(1, n_frames);
else
    if ~exist('i', 'var') || ~exist('n_steps', 'var')
        error('resume_from_workspace requires the interrupted h_atom workspace to still be loaded.');
    end
    if ~exist('active_params_set_name', 'var') || ~strcmp(params_set_name, active_params_set_name)
        error(['resume_from_workspace is true, but params_set_name does not match the active run. ' ...
               'Set resume_from_workspace=false before starting a new parameter set.']);
    end
    if i >= n_steps
        error(['resume_from_workspace is true, but the active run has already reached n_steps. ' ...
               'Set resume_from_workspace=false before starting a new parameter set.']);
    end
    if make_live_plots && (~exist('fig', 'var') || ~isgraphics(fig))
        fig = figure('Units', 'pixels', 'Position', [50, 50, 2100, 1350]);
    end
    fprintf('Resuming %s from i = %d, t = %.5e s\n', params_set_name, i + 1, (i + 1) * dt);
end

tic;
if resume_requested
    i_start = i + 1;
else
    i_start = 1;
end
for i = i_start:n_steps
    % Generate Gaussian Increments for Radial Component (r) of Brownian Motion
    dW_r = randn(M, 1) * sigma_r_factor *sigma * sqrt(dt);

    % Compute Radial Wave Function and Its Derivative at r
    R = radial_wavefunction(r, n, l);
    dR_dr = derivative(@(r) radial_wavefunction(r, n, l), r, dr);

    % Prevent Division by Zero or Very Small r Values
    r_safe = max(r, epsilon_r);

    % Handle R_safe to Prevent Division by Zero
    R_safe = max(abs(R), epsilon_R) .* (1 - 2 * (R < 0)); % sign(R);

    % Compute the physical forward drift and the geometric spherical drift
    b_r_drift = (hbar / m_e) * (dR_dr ./ R_safe);
    mu_r_curv = (hbar / m_e) * (1 ./ r_safe);
    
    % Limit the drift velocity to ensure polar energy convergence
    b_r_drift_idx = abs(b_r_drift) > v_max;
    cutoff_engagement_r_count_per_traj = cutoff_engagement_r_count_per_traj + b_r_drift_idx;
    b_r_drift(b_r_drift_idx) = v_max * sign(b_r_drift(b_r_drift_idx));

    % Limit the drift velocity (due to curviture) to ensure radial energy convergence
    mu_r_curv_idx = abs(mu_r_curv) > v_max;
    mu_r_curv(mu_r_curv_idx) = v_max * sign(mu_r_curv(mu_r_curv_idx));

    % Update the Radial Position
    r = r + (mu_r_curv + b_r_drift) .* dt + dW_r;

    % Reflective Boundary at r = 0
    r = abs(r);

    if use_radial_node_diagnostic
        current_node_value = r - 2 * a_0;
        N_cross_per_traj = N_cross_per_traj + (previous_node_value .* current_node_value < 0);
        previous_node_value = current_node_value;
        occupation_count_per_traj = occupation_count_per_traj + (r < 2 * a_0);
    end

    if make_live_plots
        % Accumulate Histogram Counts for R
        if single_traj_main_plots
            hist_counts_r = hist_counts_r + histcounts(r(1), hist_bins_r);
        else
            hist_counts_r = hist_counts_r + histcounts(r, hist_bins_r);
        end
    end

    % Generate Gaussian Increments for Brownian Motion
    dW_theta = (sigma * sqrt(dt)) .* randn(M, 1);

    % Prevent theta from Reaching Exactly 0 or pi
    theta_safe = max(min(theta, pi - epsilon_theta), epsilon_theta);

    % Compute Polar Wavefunction and Its Derivative at theta
    Y_theta = angular_wavefunction_theta(theta_safe, l, m);
    dY_theta_dtheta = derivative(@(theta) angular_wavefunction_theta(theta, l, m), theta_safe, dtheta);

    % Handle Y_theta_safe to Prevent Division by Zero
    Y_theta_safe = max(abs(Y_theta), epsilon_Y_theta) .* (1 - 2 * (Y_theta < 0)); % sign(Y_theta);

    % Compute Polar Forward Drift
    b_theta_drift = (hbar / m_e) * (dY_theta_dtheta ./ Y_theta_safe) ./ r_safe;

    % Limit the drift velocity to ensure polar energy convergence
    b_theta_drift_idx = abs(b_theta_drift) > v_max;
    cutoff_engagement_theta_count_per_traj = cutoff_engagement_theta_count_per_traj + b_theta_drift_idx;
    b_theta_drift(b_theta_drift_idx) = v_max * sign(b_theta_drift(b_theta_drift_idx));

    mu_theta_curv = (hbar / m_e) * (cot(theta_safe) / 2 ) ./ r_safe;

    % Limit the drift velocity (due to curviture) to ensure polar energy convergence
    mu_theta_curv_idx = abs(mu_theta_curv) > v_max;
    mu_theta_curv(mu_theta_curv_idx) = v_max * sign(mu_theta_curv(mu_theta_curv_idx));
   
    % Update the Polar Angle
    theta = theta + (1./ r_safe) .* ((mu_theta_curv + b_theta_drift) .* dt + dW_theta);

    if use_polar_node_diagnostic
        current_node_value = cos(theta);
        N_cross_per_traj = N_cross_per_traj + (previous_node_value .* current_node_value < 0);
        previous_node_value = current_node_value;
        occupation_count_per_traj = occupation_count_per_traj + (theta < pi / 2);
    end
  
    if make_live_plots
        % Accumulate Histogram Counts for theta
        if single_traj_main_plots
            hist_counts_theta = hist_counts_theta + histcounts(theta(1), hist_bins_theta);
        else
            hist_counts_theta = hist_counts_theta + histcounts(theta, hist_bins_theta);
        end
    end
    
    % Generate Random Steps (Gaussian Increments for Brownian Motion)
    dW_phi = (sigma * sqrt(dt)) .* randn(M, 1);

    % Compute the azimuthal forward drift from the phase Arg(Phi_m)=m*phi
    b_phi_drift = (hbar / m_e) * m ./ (r_safe .* sin(theta_safe));

    % Limit the drift velocity to ensure polar energy convergence
    b_phi_drift_idx = abs(b_phi_drift) > v_max;
    cutoff_engagement_phi_count_per_traj = cutoff_engagement_phi_count_per_traj + b_phi_drift_idx;
    b_phi_drift(b_phi_drift_idx) = v_max * sign(b_phi_drift(b_phi_drift_idx));

    % Update Azimuthal Angle
    phi_unwrapped = phi_unwrapped + 1./(r_safe .* sin(theta_safe)) .* (b_phi_drift .* dt + dW_phi);

    % Adjust phi to be Within [0, 2*pi]
    phi = mod(phi_unwrapped, 2 * pi);

    if make_live_plots
        % Accumulate Histogram Counts for phi
        if single_traj_main_plots
            hist_counts_phi = hist_counts_phi + histcounts(phi(1), hist_bins_phi);
        else
            hist_counts_phi = hist_counts_phi + histcounts(phi, hist_bins_phi);
        end
    end

    % Store only the first particle's trajectory for 3D visualisation.
    if make_live_plots && i <= max_traj
        r_traj(i)     = r(1);
        theta_traj(i) = theta(1);
        phi_traj(i)   = phi(1);
    end

    % Compute kinetic energy using the physical forward drift components
    KE_radial = 0.5 * m_e *  (b_r_drift.^2); 
    KE_theta = 0.5 * m_e * (b_theta_drift.^2); 
    KE_phi = 0.5 * m_e * (b_phi_drift.^2); 

    % Compute Potential Energy
    V = -e_charge^2 ./ (4 * pi * epsilon_0 * r_safe);
    
    % Compute Total Energy
    E = KE_radial + KE_theta + KE_phi + V;
    
    % Accumulate Energies (scalar accumulators feed the main-plot running
    % averages in subplots 6 and 7).
    if single_traj_main_plots
        S_KE_radial = S_KE_radial + KE_radial(1);
        S_KE_theta  = S_KE_theta  + KE_theta(1);
        S_KE_phi    = S_KE_phi    + KE_phi(1);
        S_b_phi     = S_b_phi     + b_phi_drift(1);
        S_V         = S_V         + V(1);
        S_E         = S_E         + E(1);
    else
        S_KE_radial = S_KE_radial + mean(KE_radial);
        S_KE_theta  = S_KE_theta  + mean(KE_theta);
        S_KE_phi    = S_KE_phi    + mean(KE_phi);
        S_b_phi     = S_b_phi     + mean(b_phi_drift);
        S_V         = S_V         + mean(V);
        S_E         = S_E         + mean(E);
    end
    % Per-particle accumulators (size [M,1]) for end-of-run SEM and for
    % the row-3 live per-trajectory diagnostic histograms. Always on,
    % independent of single_traj_main_plots.
    S_KE_radial_M = S_KE_radial_M + KE_radial;
    S_KE_theta_M  = S_KE_theta_M  + KE_theta;
    S_KE_phi_M    = S_KE_phi_M    + KE_phi;
    S_V_M         = S_V_M         + V;
    S_E_M         = S_E_M         + E;

    if mod(i, frame_step) == 0
        fprintf('i=%d, t=%.5e, KE_r=%.5e, K_azim=%.5e, V=%.5e, E=%.5e, N_cross=%d\r', ...
            i, i*dt, S_KE_radial/i, S_KE_theta / i, S_V/i, S_E/i, sum(N_cross_per_traj));

        KE_radial_traj(frame_index) = S_KE_radial /  i;
        KE_theta_traj(frame_index) = S_KE_theta /  i;
        KE_phi_traj(frame_index) = S_KE_phi / i;
        b_phi_avg_traj(frame_index) = S_b_phi / i;
        V_traj(frame_index) = S_V / i;
        E_traj(frame_index) = S_E / i;
        N_cross_total_traj(frame_index) = sum(N_cross_per_traj);
        N_cross_mean_traj(frame_index) = mean(N_cross_per_traj);
        % Per-particle running averages at this frame.
        KE_radial_traj_M(:, frame_index) = S_KE_radial_M / i;
        KE_theta_traj_M(:, frame_index)  = S_KE_theta_M  / i;
        KE_phi_traj_M(:, frame_index)    = S_KE_phi_M    / i;
        V_traj_M(:, frame_index)         = S_V_M         / i;
        E_traj_M(:, frame_index)         = S_E_M         / i;

        if n == 2 && l == 1 && m == 0
            kinetic_per_traj_now = S_KE_theta_M / i;
        elseif n == 2 && l == 0 && m == 0
            kinetic_per_traj_now = S_KE_radial_M / i;
        else
            kinetic_per_traj_now = S_E_M / i;
        end
        if ~isnan(kinetic_scan_ref) && kinetic_scan_ref ~= 0
            kinetic_mean_relerr_traj(frame_index) = abs(mean(kinetic_per_traj_now) - kinetic_scan_ref) / abs(kinetic_scan_ref);
            kinetic_median_relerr_traj(frame_index) = abs(median(kinetic_per_traj_now) - kinetic_scan_ref) / abs(kinetic_scan_ref);
            kinetic_sem_norm_traj(frame_index) = std(kinetic_per_traj_now) / sqrt(length(kinetic_per_traj_now)) / abs(kinetic_scan_ref);
            kinetic_iqr_norm_traj(frame_index) = ...
                (quantile(kinetic_per_traj_now, 0.75) - quantile(kinetic_per_traj_now, 0.25)) / abs(kinetic_scan_ref);
        end

        if make_live_plots
            hist_counts_r_traj(frame_index, :) = hist_counts_r;
            hist_counts_theta_traj(frame_index, :) = hist_counts_theta;
            hist_counts_phi_traj(frame_index, :) = hist_counts_phi;
    
            dev_arr_R(frame_index) = compute_distributions_deviation(...
                hist_bins_r, hist_counts_r, r_values, P_analytic_r);

            dev_arr_theta(frame_index) = compute_distributions_deviation(...
                hist_bins_theta, hist_counts_theta, theta_values, P_theta_analytic);

            dev_arr_phi(frame_index) = compute_distributions_deviation(...
                    hist_bins_phi, hist_counts_phi, phi_values, P_phi_analytic);        
        
            if ~isgraphics(fig)
                fig = figure('Units', 'pixels', 'Position', [50, 50, 2100, 1350]);
            end
            clf(fig);
            idx = 1:frame_index;
            if is_kinetic_scan_run && ~isnan(kinetic_scan_ref) && kinetic_scan_ref ~= 0
                subplot(2,4,1);
                plot_radial_distribution(hist_bins_r, hist_counts_r, r_values, P_analytic_r);
                subplot(2,4,2);
                plot_polar_distribution(hist_bins_theta, hist_counts_theta, theta_values, P_theta_analytic);
                subplot(2,4,3);
                plot_azimuthal_distribution(hist_bins_phi, hist_counts_phi, phi_values, P_phi_analytic);

                subplot(2,4,5);
                plot_distributions_deviation(time_arr(idx), dev_arr_R(idx), dev_arr_theta(idx), dev_arr_phi(idx));
                subplot(2,4,6);
                plot_energies_with_bands(n, l, m, time_arr(idx), ...
                    KE_radial_traj_M(:, idx), KE_theta_traj_M(:, idx), KE_phi_traj_M(:, idx), ...
                    V_traj_M(:, idx), E_traj_M(:, idx), false);

                subplot(2,4,7);
                plot_azimuthal_current(time_arr(idx), b_phi_avg_traj(idx), n, l, m);

                subplot(2,4,4);
                hold on;
                plot(time_arr(idx), kinetic_mean_relerr_traj(idx), ...
                    'Color', [0.00, 0.25, 0.55], 'LineWidth', 1.2);
                plot(time_arr(idx), kinetic_median_relerr_traj(idx), ...
                    'Color', [0.55, 0.25, 0.00], 'LineWidth', 1.2);
                hold off;
                xlabel('Time (s)');
                ylabel('relative error', 'Interpreter', 'latex');
                title(['Accuracy of ' kinetic_scan_label], 'Interpreter', 'latex');
                legend({'mean', 'median'}, 'Interpreter', 'latex', ...
                       'Location', 'northeast', 'Box', 'off');
                grid on;

                subplot(2,4,8);
                hold on;
                plot(time_arr(idx), kinetic_sem_norm_traj(idx), ...
                    'Color', [0.45, 0.10, 0.55], 'LineWidth', 1.2);
                plot(time_arr(idx), kinetic_iqr_norm_traj(idx), ...
                    'Color', [0.10, 0.45, 0.25], 'LineWidth', 1.2);
                hold off;
                xlabel('Time (s)');
                ylabel('normalised spread', 'Interpreter', 'latex');
                title(['Trajectory spread of ' kinetic_scan_label], 'Interpreter', 'latex');
                legend({'SEM$/T_{\rm an}$', 'IQR$/T_{\rm an}$'}, 'Interpreter', 'latex', ...
                       'Location', 'northeast', 'Box', 'off');
                grid on;
            else
                subplot(3,4,1);
                plot_radial_distribution(hist_bins_r, hist_counts_r, r_values, P_analytic_r);
                subplot(3,4,2);
                plot_polar_distribution(hist_bins_theta, hist_counts_theta, theta_values, P_theta_analytic);
                subplot(3,4,3);
                plot_azimuthal_distribution(hist_bins_phi, hist_counts_phi, phi_values, P_phi_analytic);
                subplot(3,4,4);
                plot_trajectory_3d(n, l, r_traj, theta_traj, phi_traj, a_0, i, params_set_name);

                subplot(3,4,5);
                plot_distributions_deviation(time_arr(idx), dev_arr_R(idx), dev_arr_theta(idx), dev_arr_phi(idx));
                subplot(3,4,6);
                if single_traj_main_plots
                    plot_energies(n, l, m, time_arr(idx), KE_radial_traj(idx), ...
                        KE_theta_traj(idx), KE_phi_traj(idx), V_traj(idx), E_traj(idx));
                else
                    plot_energies_with_bands(n, l, m, time_arr(idx), ...
                        KE_radial_traj_M(:, idx), KE_theta_traj_M(:, idx), KE_phi_traj_M(:, idx), ...
                        V_traj_M(:, idx), E_traj_M(:, idx));
                end

                subplot(3,4,7);
                plot_azimuthal_current(time_arr(idx), b_phi_avg_traj(idx), n, l, m);

                subplot(3,4,8);
                plot_trajectory_xy(n, l, r_traj, theta_traj, phi_traj, a_0, i, params_set_name);

                % Row 3: per-trajectory energy summaries across the M
                % independent trajectories.
                subplot(3,4,9);
                plot_per_traj_energy_summary_panel(S_KE_radial_M / i, ...
                    '$\left< E_r \right>$', ref_energies.E_r);
                subplot(3,4,10);
                plot_per_traj_energy_summary_panel(S_KE_theta_M / i, ...
                    '$\left< E_\theta \right>$', ref_energies.E_theta);
                subplot(3,4,11);
                plot_per_traj_energy_summary_panel(S_V_M / i, ...
                    '$\left< V \right>$', ref_energies.V);
                subplot(3,4,12);
                plot_per_traj_energy_summary_panel(S_E_M / i, ...
                    '$\left< E \right>$', ref_energies.E);
            end

            annotation_text = sprintf('t = %.2e s', time_arr(frame_index));
            annotation('textbox', [0.02, 0.95, 0.1, 0.05], 'String', annotation_text, ...
            'Units', 'normalized', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', ...
            'FontSize', 12, 'FontWeight', 'bold', 'EdgeColor', 'none');

            drawnow;
            if isgraphics(fig)
                frame = getframe(fig);
                % VideoWriter locks the frame size on the first frame. On
                % HiDPI/Retina displays getframe can return a slightly
                % different pixel size mid-run (window moved / DPI rescaled),
                % which would abort the whole run. Rescale to the locked
                % size so the simulation never dies on a video write.
                if ~isempty(v.Height) && ...
                        (size(frame.cdata, 1) ~= v.Height || size(frame.cdata, 2) ~= v.Width)
                    frame.cdata = imresize(frame.cdata, [v.Height, v.Width]);
                end
                writeVideo(v, frame);
            end
        end
        frame_index = frame_index + 1;
    end
end

totalElapsedTime = toc;
if make_live_plots
    close(v);
end
fprintf('Total elapsed time: %.4f seconds\n', totalElapsedTime);

if make_live_plots && isgraphics(fig)
    exportgraphics(fig, fullfile(base_dir, strcat(output_run_name,'.pdf')), 'ContentType', 'image');
    clear fig;
end

% =====================================================================
% Per-trajectory time-averaged energies and statistical uncertainty
% (see Section 11 of the manuscript).
% =====================================================================
KE_r_per_traj = S_KE_radial_M / n_steps;
KE_t_per_traj = S_KE_theta_M  / n_steps;
KE_p_per_traj = S_KE_phi_M    / n_steps;
V_per_traj    = S_V_M         / n_steps;
E_per_traj    = S_E_M         / n_steps;

f_occupation_per_traj = occupation_count_per_traj / n_steps;
cutoff_engagement_r_per_traj = cutoff_engagement_r_count_per_traj / n_steps;
cutoff_engagement_theta_per_traj = cutoff_engagement_theta_count_per_traj / n_steps;
cutoff_engagement_phi_per_traj = cutoff_engagement_phi_count_per_traj / n_steps;
N_cross_total = sum(N_cross_per_traj);
N_cross_mean = mean(N_cross_per_traj);
N_cross_std = std(N_cross_per_traj);

report_sem = @(label, x) fprintf('  %-12s = %+.4e +/- %.2e J   (M=%d traj)\n', ...
    label, mean(x), std(x)/sqrt(length(x)), length(x));
report_med = @(label, x) fprintf('  %-12s    median = %+.4e ,  IQR = [%+.4e, %+.4e] J\n', ...
    label, median(x), quantile(x, 0.25), quantile(x, 0.75));

fprintf('\n=== Per-trajectory energy averages (state %s) ===\n', params_set_name);
report_sem('<E_r>',       KE_r_per_traj);
report_sem('<E_theta>',   KE_t_per_traj);
report_sem('<E_phi>',     KE_p_per_traj);
report_sem('<V>',         V_per_traj);
report_sem('<E>',         E_per_traj);
fprintf('--- Robust statistics for the (potentially) heavy-tailed observables ---\n');
report_med('<E_theta>',   KE_t_per_traj);
report_med('<E>',         E_per_traj);
fprintf('============================================\n');

fprintf('\n=== Cutoff and nodal diagnostics (state %s, v_max = %.4g c) ===\n', ...
        params_set_name, v_max_over_c);
if use_polar_node_diagnostic
    fprintf('  Nodal observable: cos(theta) sign changes; occupation = fraction(theta < pi/2)\n');
elseif use_radial_node_diagnostic
    fprintf('  Nodal observable: r - 2*a_0 sign changes; occupation = fraction(r < 2*a_0)\n');
else
    fprintf('  No nodal occupation diagnostic is defined for this state.\n');
end
fprintf('  %-5s %10s %14s %12s %12s %12s\n', ...
        'traj', 'N_cross', 'f_occupation', 'cutoff_r', 'cutoff_th', 'cutoff_phi');
for j = 1:M
    fprintf('  %-5d %10.0f %14.5f %12.3e %12.3e %12.3e\n', ...
            j, N_cross_per_traj(j), f_occupation_per_traj(j), ...
            cutoff_engagement_r_per_traj(j), cutoff_engagement_theta_per_traj(j), ...
            cutoff_engagement_phi_per_traj(j));
end
fprintf('  total crossings = %.0f, mean/traj = %.2f, std/traj = %.2f\n', ...
        N_cross_total, N_cross_mean, N_cross_std);
fprintf('===============================================================\n');

resume_from_workspace = false;
resume_requested = false;
save(fullfile(base_dir, strcat(output_run_name,'.mat')), '-v7.3');

function plot_per_traj_energy_summary_panel(values, label, ref_value)
    mean_value = mean(values);
    sem_value = std(values) / sqrt(length(values));
    median_value = median(values);
    iqr_value = quantile(values, 0.75) - quantile(values, 0.25);

    axis off;
    title(label, 'Interpreter', 'latex');
    text(0.05, 0.78, sprintf('mean = %.4e J', mean_value), 'Units', 'normalized');
    text(0.05, 0.60, sprintf('SEM = %.2e J', sem_value), 'Units', 'normalized');
    text(0.05, 0.42, sprintf('median = %.4e J', median_value), 'Units', 'normalized');
    text(0.05, 0.24, sprintf('IQR = %.2e J', iqr_value), 'Units', 'normalized');
    if ~isnan(ref_value)
        rel_error = abs(mean_value - ref_value) / abs(ref_value);
        text(0.05, 0.06, sprintf('rel.err = %.3g', rel_error), 'Units', 'normalized');
    end
end
