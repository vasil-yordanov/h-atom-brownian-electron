function make_dashboard_movie(run_name, fig_width, frame_stride)
% MAKE_DASHBOARD_MOVIE  Rebuild the 6-panel dashboard movie from saved data,
% WITHOUT re-running the simulation.
%
% This is a replica of the live per-frame rendering loop in h_atom.m: it reads
% the per-frame arrays that the simulation already stored in
% data/<run_name>/<run_name>.mat and redraws the same dashboard frame by frame.
% Use it to restyle the panels (layout, 3D box, colours, ...) and regenerate the
% movie in minutes instead of re-running the (hours-long) simulation.
%
% Usage:
%   make_dashboard_movie('1s0_1fs')
%   make_dashboard_movie('2p0_10ps')
%
% Writes data/<run_name>/<run_name>.mp4 (same path/name as the live run), so the
% existing MP4 -> GIF ffmpeg step in the README works unchanged.
%
% Notes:
%   * Only the production (6-panel, 2x3) dashboard is reproduced -- the same one
%     the state movies use. (Scan runs use a different live dashboard.)
%   * The frame content matches the live run exactly: panel histograms come from
%     hist_counts_*_traj(k,:), the time-series panels from the (1:k) slices, and
%     the 3D trajectory from traj_points = k*frame_step (column index = step
%     index, exactly as h_atom stores r_traj).

    addpath('functions');
    run_name = char(run_name);   % accept "double-quoted" strings too

    % Speed knobs. Defaults are tuned to the README GIF target (1200 px wide at
    % 5 fps), which is ~6x faster than the full live-sim replica and produces the
    % same final GIF. Pass make_dashboard_movie('1s0_1fs', 2100, 1) for the exact
    % full-fidelity replica of the simulation's own MP4.
    %   fig_width    : figure width in px (height scales to keep 2100x1350 ratio)
    %   frame_stride : render every Nth frame; the frame rate is divided by N so
    %                  the movie keeps the same duration.
    if nargin < 2 || isempty(fig_width),    fig_width = 1200; end
    if nargin < 3 || isempty(frame_stride), frame_stride = 2; end

    mat_path = fullfile('data', run_name, [run_name '.mat']);
    if ~exist(mat_path, 'file')
        error('Missing %s -- run the simulation (h_atom) for this set first.', mat_path);
    end
    S = load(mat_path);

    n = S.n;  l = S.l;  m = S.m;  a_0 = S.a_0;
    n_frames   = double(S.n_frames);
    frame_step = double(S.frame_step);
    params_set_name = S.params_set_name;        % char row, used by plot_trajectory_3d
    single_traj_main_plots = logical(S.single_traj_main_plots);

    if ~isfield(S, 'hist_counts_r_traj') || all(S.hist_counts_r_traj(:) == 0)
        error(['%s has no per-frame dashboard data (was it run with ', ...
               'make_live_plots = false?). Cannot rebuild the movie.'], mat_path);
    end

    % Video writer: same frame rate convention as the live run.
    base_fps = 10;
    if any(strcmp(params_set_name, ...
            {'2p_m1_1ps', '2p1_1ps', '2p_m0_1ps', '2p0_1ps', '2p_mn1_1ps'}))
        base_fps = 24;
    end
    out_mp4 = fullfile('data', run_name, [run_name '.mp4']);
    vw = VideoWriter(out_mp4, 'MPEG-4');
    vw.FrameRate = max(1, base_fps / frame_stride);   % keep duration with fewer frames
    open(vw);

    fig_height = round(fig_width * 1350 / 2100);
    fig = figure('Units', 'pixels', 'Position', [50, 50, fig_width, fig_height], 'Color', 'w');

    for frame_index = 1:frame_stride:n_frames
        clf(fig);
        idx = 1:frame_index;

        % --- 6-panel (2x3) production dashboard, replicated from h_atom.m ---
        subplot(2,3,1);
        plot_radial_distribution(S.hist_bins_r, S.hist_counts_r_traj(frame_index, :), ...
                                 S.r_values, S.P_analytic_r);
        subplot(2,3,2);
        plot_polar_distribution(S.hist_bins_theta, S.hist_counts_theta_traj(frame_index, :), ...
                                S.theta_values, S.P_theta_analytic);
        subplot(2,3,3);
        plot_azimuthal_distribution(S.hist_bins_phi, S.hist_counts_phi_traj(frame_index, :), ...
                                    S.phi_values, S.P_phi_analytic);

        subplot(2,3,4);
        plot_distributions_deviation(S.time_arr(idx), S.dev_arr_R(idx), ...
                                     S.dev_arr_theta(idx), S.dev_arr_phi(idx));
        subplot(2,3,5);
        if single_traj_main_plots
            plot_energies(n, l, m, S.time_arr(idx), S.KE_radial_traj(idx), ...
                S.KE_theta_traj(idx), S.KE_phi_traj(idx), S.V_traj(idx), S.E_traj(idx));
        else
            plot_energies_with_bands(n, l, m, S.time_arr(idx), ...
                S.KE_radial_traj_M(:, idx), S.KE_theta_traj_M(:, idx), ...
                S.KE_phi_traj_M(:, idx), S.V_traj_M(:, idx), S.E_traj_M(:, idx));
        end
        subplot(2,3,6);
        plot_trajectory_3d(n, l, S.r_traj, S.theta_traj, S.phi_traj, a_0, ...
                           frame_index * frame_step, params_set_name);

        annotation_text = sprintf('t = %.2e s', S.time_arr(frame_index));
        annotation('textbox', [0.02, 0.95, 0.1, 0.05], 'String', annotation_text, ...
            'Units', 'normalized', 'HorizontalAlignment', 'left', ...
            'VerticalAlignment', 'top', 'FontSize', 12, 'FontWeight', 'bold', ...
            'EdgeColor', 'none');

        drawnow;
        frame = getframe(fig);
        % Match the locked video size (HiDPI getframe can drift; same guard as h_atom).
        if ~isempty(vw.Height) && ...
                (size(frame.cdata, 1) ~= vw.Height || size(frame.cdata, 2) ~= vw.Width)
            frame.cdata = imresize(frame.cdata, [vw.Height, vw.Width]);
        end
        writeVideo(vw, frame);

        if mod(frame_index, 100) == 0
            fprintf('  %d / %d frames\n', frame_index, n_frames);
        end
    end

    close(vw);
    if isgraphics(fig); close(fig); end
    fprintf('Dashboard movie rebuilt: %s (%d frames)\n', out_mp4, n_frames);
end
