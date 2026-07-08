function make_trajectory_video(state)
% MAKE_TRAJECTORY_VIDEO  Build the final 3-D electron trajectory video.
%
% Usage:
%   make_trajectory_video("1s0")
%   make_trajectory_video("2s0")
%   make_trajectory_video("2p0")
%
% Outputs:
%   1s0 -> videos/1s0_native1200.mp4  (0.15 ps, box +/-6 a0)
%   2s0 -> videos/2s0_native1200.mp4  (0.5 ps,  box +/-16 a0)
%   2p0 -> videos/2p0_native1200.mp4  (2 ps,    box +/-16 a0)

    if nargin < 1 || isempty(state)
        error('Usage: make_trajectory_video("1s0"), make_trajectory_video("2s0"), or make_trajectory_video("2p0")');
    end

    cfg = trajectory_video_config(state);
    if ~exist(cfg.mat_path, 'file')
        error('Missing %s. Run the corresponding simulation first.', cfg.mat_path);
    end

    render_video(cfg);
end

function cfg = trajectory_video_config(state)
    key = lower(strtrim(char(state)));
    key = strrep(key, '-', '_');

    cfg.n_frames = 1000;
    cfg.az_span = 60;
    cfg.npts = 3e6;
    cfg.minlen = 8;
    cfg.crop_pct = 97;
    cfg.show_labels = true;
    cfg.cut_mode = 'verticalquarter';
    cfg.final_az = 165;
    cfg.final_el = 12;
    cfg.fig_px = 600; % HiDPI capture gives about 1200 px video width.

    switch key
        case {'1s0', '1s0_movie'}
            cfg.mat_path = 'data/1s0_movie/1s0_movie.mat';
            cfg.out = 'videos/1s0_native1200.mp4';
            cfg.box_limit = 6;

        case {'2s0', '2s0_movie'}
            cfg.mat_path = 'data/2s0_movie/2s0_movie.mat';
            cfg.out = 'videos/2s0_native1200.mp4';
            cfg.box_limit = 16;

        case {'2p0', '2p0_movie'}
            cfg.mat_path = 'data/2p0_movie/2p0_movie.mat';
            cfg.out = 'videos/2p0_native1200.mp4';
            cfg.box_limit = 16;

        otherwise
            error('Unknown state "%s". Use "1s0", "2s0", or "2p0".', state);
    end
end

function render_video(cfg)
    addpath('functions');

    mf = matfile(cfg.mat_path);
    a_0 = 5.29177e-11;
    try
        a_0 = double(mf.a_0);
    catch
    end
    dt = double(mf.dt);

    rau = single(mf.r_traj) ./ single(a_0);
    th = single(mf.theta_traj);
    ph = single(mf.phi_traj);
    x = rau .* sin(th) .* cos(ph);
    y = rau .* sin(th) .* sin(ph);
    z = rau .* cos(th);
    clear th ph

    xyz_max = double([max(abs(x(:))), max(abs(y(:))), max(abs(z(:)))]);
    color_lim = double(ceil(prctile(rau, cfg.crop_pct)));

    N = numel(x);
    s = max(1, round(N / cfg.npts));
    X = x(1:s:end);
    Y = y(1:s:end);
    Z = z(1:s:end);
    C = rau(1:s:end);
    clear x y z rau

    valid = cut_mask(X, Y, Z, cfg.cut_mode);
    keep = drop_short_runs(valid, cfg.minlen);
    X(~keep) = NaN;
    Y(~keep) = NaN;
    Z(~keep) = NaN;
    C(~keep) = NaN;
    npts_a = numel(X);

    fig_h = round(cfg.fig_px * 740 / 760);
    fig = figure('Visible', 'off', 'Color', 'w', 'Position', [50 50 cfg.fig_px fig_h]);
    ax = axes(fig);
    h = surface(ax, [X; X], [Y; Y], [Z; Z], [C; C], ...
                'FaceColor', 'none', 'EdgeColor', 'interp', ...
                'LineWidth', 0.35, 'EdgeAlpha', 0.3);
    xlim(ax, [-cfg.box_limit cfg.box_limit]);
    ylim(ax, [-cfg.box_limit cfg.box_limit]);
    zlim(ax, [-cfg.box_limit cfg.box_limit]);
    set(ax, 'DataAspectRatio', [1 1 1], 'PlotBoxAspectRatio', [1 1 1]);
    camup(ax, [0 0 1]);
    grid(ax, 'on');
    box(ax, 'on');
    colormap(ax, turbo);
    clim(ax, [0 color_lim]);
    cb = colorbar(ax);
    cb.Label.String = 'r / a_0';

    if cfg.show_labels
        xlabel(ax, 'X / a_0');
        ylabel(ax, 'Y / a_0');
        zlabel(ax, 'Z / a_0');
        time_txt = annotation(fig, 'textbox', [0.10 0.91 0.34 0.05], ...
                              'String', '', 'EdgeColor', 'none', ...
                              'Color', [0.10 0.10 0.10], 'FontSize', 11, ...
                              'FontWeight', 'bold', 'Interpreter', 'tex');
    else
        xlabel(ax, '');
        ylabel(ax, '');
        zlabel(ax, '');
        time_txt = [];
    end

    out_dir = fileparts(cfg.out);
    if ~isempty(out_dir) && ~exist(out_dir, 'dir')
        mkdir(out_dir);
    end
    vw = VideoWriter(cfg.out, 'MPEG-4');
    vw.FrameRate = 30;
    open(vw);

    for k = 1:cfg.n_frames
        m = max(2, round(k / cfg.n_frames * npts_a));
        Xk = X;
        Xk(m+1:end) = NaN;
        set(h, 'XData', [Xk; Xk]);

        az = cfg.final_az - cfg.az_span + (k - 1) / (cfg.n_frames - 1) * cfg.az_span;
        view(ax, az, cfg.final_el);
        if cfg.show_labels
            set(time_txt, 'String', trajectory_time_label(m * s * dt));
        end
        drawnow;
        writeVideo(vw, getframe(fig));
        if mod(k, 100) == 0
            fprintf('  %d / %d\n', k, cfg.n_frames);
        end
    end

    close(vw);
    close(fig);

    excess = max(0, xyz_max - cfg.box_limit);
    if any(excess > 0)
        clip_desc = sprintf(', clipping beyond box [%g %g %g] a0', excess(1), excess(2), excess(3));
    else
        clip_desc = ', no coordinate clipping';
    end
    fprintf('wrote %s (%d frames, az %g->%g, el %g, cut %s, minlen %d, color p%g, xyz max [%g %g %g], box fixed +/- %g a0%s)\n', ...
            cfg.out, cfg.n_frames, cfg.final_az - cfg.az_span, cfg.final_az, cfg.final_el, ...
            cfg.cut_mode, cfg.minlen, cfg.crop_pct, xyz_max(1), xyz_max(2), xyz_max(3), ...
            cfg.box_limit, clip_desc);
end

function label = trajectory_time_label(seconds)
    if ~isfinite(seconds)
        label = 't = n/a';
        return;
    end
    fs = seconds / 1e-15;
    label = sprintf('t = %.0f fs', fs);
end

function valid = cut_mask(X, Y, Z, cut_mode)
    switch lower(char(cut_mode))
        case 'octant'
            valid = ~(X > 0 & Y > 0 & Z > 0);
        case {'verticalquarter', 'vertical-quarter', 'vertical_quarter', 'quarter', 'quarterspace', 'quarter-space'}
            valid = ~(X > 0 & Y > 0);
        case {'vertical', 'xhalfspace', 'x-halfspace', 'x_halfspace'}
            valid = ~(X > 0);
        case {'halfspace', 'half-space', 'half_space', 'yhalfspace', 'y-halfspace', 'y_halfspace'}
            valid = ~(Y > 0);
        case {'quadrupole', 'quadropol', 'quad'}
            valid = ~(X .* Y .* Z > 0);
        case {'none', 'full', 'nocut', 'no-cut', 'no_cut'}
            valid = true(size(X));
        otherwise
            error('Unknown cut_mode "%s". Use octant, verticalquarter, vertical, halfspace, quadrupole, or none.', cut_mode);
    end
end

function keep = drop_short_runs(valid, minlen)
    if minlen <= 1
        keep = logical(valid(:)');
        return;
    end
    v = double(valid(:)');
    starts = find(diff([0 v]) == 1);
    ends = find(diff([v 0]) == -1);
    keep = logical(v);
    for i = find((ends - starts + 1) < minlen)
        keep(starts(i):ends(i)) = false;
    end
end
