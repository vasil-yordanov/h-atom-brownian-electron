function make_dashboard_movie(run_name, fig_width, frame_stride, max_time_override)
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
%   make_dashboard_movie('sphere')
%   make_dashboard_movie('1s0')
%   make_dashboard_movie('2p0')
%   make_dashboard_movie('1s0_movie', [], [], 0.05e-12)   % cut at 0.05 ps
%
% Pass max_time_override (4th arg, seconds) to trim where the movie stops without
% re-running the simulation -- it just re-encodes fewer frames. Use it to pick the
% cut time off the rendered video.
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
    requested_run_name = char(run_name);   % accept "double-quoted" strings too
    trajectory_style = 'radius';
    switch requested_run_name
        case 'sphere'
            run_name = '1s0_1fs';
            trajectory_style = 'sphere';
        case '1s0'
            run_name = '1s0_1ps';
        case '2s0'
            run_name = '2s0_10ps';
        case '2p0'
            run_name = '2p0_10ps';
        case '2p_m1'
            run_name = '2p_m1_1ps';
        otherwise
            run_name = requested_run_name;
    end

    % Defaults: every frame (frame_stride = 1) at the README GIF width (1200 px,
    % downsized from the 2100 px live MP4). Pass a larger frame_stride to render
    % faster -- every Nth frame, with the frame rate divided by N so the duration is
    % unchanged. Use make_dashboard_movie('1s0_1fs', 2100, 1) for the full-resolution
    % exact replica of the simulation's own MP4.
    %   fig_width    : figure width in px (height scales to keep 2100x1350 ratio)
    %   frame_stride : render every Nth frame (default 1 = keep all frames)
    if nargin < 2 || isempty(fig_width)
        fig_width = 1200;   % standard README dashboard size (downsized from the 2100 px live MP4)
    end
    if nargin < 3 || isempty(frame_stride), frame_stride = 1; end
    if nargin < 4, max_time_override = []; end   % seconds; overrides the per-state cut

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

    % Per-state time window: render the movie only up to this wall-clock time and
    % stop there (the later motion is slow / less informative). Matched by the
    % quantum numbers (n,l,m) of the run, so it applies to every run-name alias of
    % a state. inf = render the full run.
    max_time = inf;
    if     n == 1 && l == 0 && m == 0,       max_time = 0.2e-12;   % (1,0,0): 0.2 ps
    elseif n == 2 && l == 1 && m == 0,       max_time = 2.0e-12;   % (2,1,0): 2 ps
    elseif n == 2 && l == 1 && abs(m) == 1,  max_time = 2.0e-12;   % (2,1,+-1): 2 ps (1 ps run -> full)
    elseif n == 2 && l == 0 && m == 0,       max_time = 1.0e-12;   % (2,0,0): 1 ps
    end
    if ~isempty(max_time_override), max_time = max_time_override; end
    n_frames_movie = n_frames;
    if isfinite(max_time)
        last = find(S.time_arr(1:n_frames) <= max_time, 1, 'last');
        if ~isempty(last), n_frames_movie = last; end
    end

    % Video writer: same frame rate convention as the live run.
    base_fps = 10;
    if any(strcmp(params_set_name, ...
            {'2p_m1_1ps', '2p1_1ps', '2p_m0_1ps', '2p0_1ps', '2p_mn1_1ps'}))
        base_fps = 24;
    elseif any(strcmp(params_set_name, ...
            {'1s0_movie', '2p0_movie', '2s0_movie', '2p_m1_movie'}))
        base_fps = 50;   % 1000-frame movie runs -> ~20 s at 50 fps
    end
    out_mp4 = fullfile('data', run_name, [run_name '.mp4']);
    vw = VideoWriter(out_mp4, 'MPEG-4');
    vw.FrameRate = max(1, base_fps / frame_stride);   % keep duration with fewer frames
    open(vw);

    fig_height = round(fig_width * 1350 / 2100);
    fig = figure('Units', 'pixels', 'Position', [50, 50, fig_width, fig_height], 'Color', 'w');

    for frame_index = 1:frame_stride:n_frames_movie
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
                           frame_index * frame_step, params_set_name, false, ...
                           trajectory_style);

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
            fprintf('  %d / %d frames\n', frame_index, n_frames_movie);
        end
    end

    close(vw);
    if isgraphics(fig); close(fig); end
    fprintf('Dashboard movie rebuilt: %s (%d of %d frames, up to t = %.3g s)\n', ...
            out_mp4, n_frames_movie, n_frames, S.time_arr(n_frames_movie));
end
