function make_manuscript_trajectory_figures(kind, t_final)
% MAKE_MANUSCRIPT_TRAJECTORY_FIGURES  Build manuscript 3-D trajectory PDFs.
%
% Usage:
%   make_manuscript_trajectory_figures("sphere", 200e-18)
%   make_manuscript_trajectory_figures("sphere", 500e-18)
%   make_manuscript_trajectory_figures("2s0")
%   make_manuscript_trajectory_figures("2p0")
%
% For "sphere", numeric t_final values smaller than 1e-12 are interpreted as
% seconds; larger numeric values are interpreted as attoseconds. Strings such as
% "200as", "0.2fs", and "2e-16s" are also accepted.

    if nargin < 1 || isempty(kind)
        error(['Usage: make_manuscript_trajectory_figures("sphere", t_final), ', ...
               'make_manuscript_trajectory_figures("2s0"), or ', ...
               'make_manuscript_trajectory_figures("2p0")']);
    end

    key = lower(strtrim(char(kind)));
    key = strrep(key, '-', '_');

    switch key
        case {'sphere', '1s0_1fs', 'fig1'}
            if nargin < 2 || isempty(t_final)
                error('Usage: make_manuscript_trajectory_figures("sphere", t_final)');
            end
            cfg = sphere_config();
            times = parse_time_values(t_final);
            for seconds = times
                render_sphere_panel(cfg, seconds);
            end

        case {'1s0', '1s0_movie', '2s0', '2s0_movie', '2p0', '2p0_movie'}
            if nargin > 1 && ~isempty(t_final)
                error('The "%s" panel uses its stored final time; omit t_final.', kind);
            end
            cfg = cutaway_config(key);
            render_cutaway_panel(cfg);

        otherwise
            error('Unknown panel "%s". Use "sphere", "1s0", "2s0", or "2p0".', kind);
    end
end

function cfg = sphere_config()
    cfg.mat_path = 'data/1s0_1fs/1s0_1fs.mat';
    cfg.out_dir = 'figures';
    cfg.run_name = '1s0_1fs';
    cfg.view_az = -25;
    cfg.view_el = 10;
    cfg.export_resolution = 600;
    cfg.font_scale = 2.5;
    cfg.line_color = [0.10 0.35 0.80];
    cfg.line_width = 0.9;
    cfg.line_alpha = 0.75;
    cfg.ax_position = [0.16 0.18 0.76 0.70];
    cfg.time_position = [0.16 0.85 0.34 0.06];
end

function cfg = cutaway_config(key)
    cfg.npts = 3e6;
    cfg.minlen = 8;
    cfg.crop_pct = 97;
    cfg.font_scale = 2.5;
    cfg.cut_mode = 'verticalquarter';
    cfg.view_az = 165;
    cfg.view_el = 12;
    cfg.export_resolution = 600;
    cfg.ax_position = [0.15 0.18 0.61 0.70];
    cfg.colorbar_position = [0.86 0.20 0.025 0.64];
    cfg.time_position = [0.15 0.85 0.34 0.06];

    switch key
        case {'1s0', '1s0_movie'}
            cfg.mat_path = 'data/1s0_movie/1s0_movie.mat';
            cfg.out_base = 'figures/1s0_0p15ps';
            cfg.box_limit = 6;

        case {'2s0', '2s0_movie'}
            cfg.mat_path = 'data/2s0_movie/2s0_movie.mat';
            cfg.out_base = 'figures/2s0_0p5ps';
            cfg.box_limit = 16;

        case {'2p0', '2p0_movie'}
            cfg.mat_path = 'data/2p0_movie/2p0_movie.mat';
            cfg.out_base = 'figures/2p0_2ps';
            cfg.box_limit = 16;
    end
end

function render_sphere_panel(cfg, requested_seconds)
    addpath('functions');
    require_file(cfg.mat_path, 'Run the 1s0_1fs simulation first.');

    mf = matfile(cfg.mat_path);
    n = double(mf.n);
    l = double(mf.l);
    dt = double(mf.dt);
    a_0 = load_a0(mf);

    sz = size(mf, 'r_traj');
    N = sz(2);
    npts = min(max(2, round(requested_seconds / dt)), N);
    actual_seconds = npts * dt;

    [X, Y, Z] = load_cartesian_prefix(mf, a_0, npts);

    fig = figure('Visible', 'off', 'Color', 'w', 'Position', [50 50 900 860]);
    ax = axes(fig);
    surface(ax, [X; X], [Y; Y], [Z; Z], ...
            'FaceColor', 'none', 'EdgeColor', cfg.line_color, ...
            'LineWidth', cfg.line_width, 'EdgeAlpha', cfg.line_alpha);
    clear X Y Z

    lim = cloud_cube_limit(n, l, a_0, cfg.run_name);
    finish_axes(ax, cfg, [-lim lim], [-lim lim], [-lim lim]);
    add_time_label(fig, cfg.time_position, time_label(actual_seconds, 'as'), 24);

    out_path = fullfile(cfg.out_dir, sphere_filename(actual_seconds));
    write_pdf(fig, out_path, cfg.export_resolution);
    fprintf('wrote %s (%d points)\n', out_path, npts);
end

function render_cutaway_panel(cfg)
    addpath('functions');
    require_file(cfg.mat_path, 'Run the corresponding simulation first.');

    mf = matfile(cfg.mat_path);
    dt = double(mf.dt);
    a_0 = load_a0(mf);

    rau = single(mf.r_traj) ./ single(a_0);
    th = single(mf.theta_traj);
    ph = single(mf.phi_traj);
    Xfull = rau .* sin(th) .* cos(ph);
    Yfull = rau .* sin(th) .* sin(ph);
    Zfull = rau .* cos(th);
    clear th ph

    xyz_max = double([max(abs(Xfull(:))), max(abs(Yfull(:))), max(abs(Zfull(:)))]);
    color_lim = double(ceil(prctile(rau, cfg.crop_pct)));

    N = numel(Xfull);
    s = max(1, round(N / cfg.npts));
    X = Xfull(1:s:end);
    Y = Yfull(1:s:end);
    Z = Zfull(1:s:end);
    C = rau(1:s:end);
    clear Xfull Yfull Zfull rau

    valid = cut_mask(X, Y, Z, cfg.cut_mode);
    keep = drop_short_runs(valid, cfg.minlen);
    X(~keep) = NaN;
    Y(~keep) = NaN;
    Z(~keep) = NaN;
    C(~keep) = NaN;

    fig = figure('Visible', 'off', 'Color', 'w', 'Position', [50 50 900 860]);
    ax = axes(fig);
    surface(ax, [X; X], [Y; Y], [Z; Z], [C; C], ...
            'FaceColor', 'none', 'EdgeColor', 'interp', ...
            'LineWidth', 0.35, 'EdgeAlpha', 0.3);
    clear X Y Z C

    finish_axes(ax, cfg, [-cfg.box_limit cfg.box_limit], ...
                [-cfg.box_limit cfg.box_limit], [-cfg.box_limit cfg.box_limit]);
    colormap(ax, turbo);
    clim(ax, [0 color_lim]);
    cb = colorbar(ax);
    cb.Position = cfg.colorbar_position;
    cb.FontSize = ax.FontSize;
    cb.Label.String = 'r / a_0';
    cb.Label.FontSize = ax.XLabel.FontSize;
    add_time_label(fig, cfg.time_position, trajectory_time_label(N * dt), 24);

    write_pdf(fig, [cfg.out_base '.pdf'], cfg.export_resolution);

    excess = max(0, xyz_max - cfg.box_limit);
    if any(excess > 0)
        clip_desc = sprintf(', clipping beyond box [%g %g %g] a0', excess(1), excess(2), excess(3));
    else
        clip_desc = ', no coordinate clipping';
    end
    fprintf('wrote %s.pdf (npts=%g, color p%g, xyz max [%g %g %g], box fixed +/- %g a0%s)\n', ...
            cfg.out_base, cfg.npts, cfg.crop_pct, xyz_max(1), xyz_max(2), xyz_max(3), ...
            cfg.box_limit, clip_desc);
end

function finish_axes(ax, cfg, x_limits, y_limits, z_limits)
    xlim(ax, x_limits);
    ylim(ax, y_limits);
    zlim(ax, z_limits);
    set(ax, 'DataAspectRatio', [1 1 1], 'PlotBoxAspectRatio', [1 1 1]);
    ax_base_font = ax.FontSize;
    camup(ax, [0 0 1]);
    grid(ax, 'on');
    box(ax, 'on');
    set(ax, 'Position', cfg.ax_position);
    view(ax, cfg.view_az, cfg.view_el);
    ax.FontSize = ax_base_font * cfg.font_scale;
    xlabel(ax, 'X / a_0');
    ylabel(ax, 'Y / a_0');
    zlabel(ax, 'Z / a_0');
    ax.XLabel.FontSize = ax_base_font * cfg.font_scale;
    ax.YLabel.FontSize = ax_base_font * cfg.font_scale;
    ax.ZLabel.FontSize = ax_base_font * cfg.font_scale;
    title(ax, '');
end

function [X, Y, Z] = load_cartesian_prefix(mf, a_0, npts)
    r = single(mf.r_traj(1, 1:npts));
    th = single(mf.theta_traj(1, 1:npts));
    ph = single(mf.phi_traj(1, 1:npts));
    rau = r ./ single(a_0);
    X = rau .* sin(th) .* cos(ph);
    Y = rau .* sin(th) .* sin(ph);
    Z = rau .* cos(th);
end

function a_0 = load_a0(mf)
    a_0 = 5.29177e-11;
    try
        a_0 = double(mf.a_0);
    catch
    end
end

function add_time_label(fig, position, label, font_size)
    annotation(fig, 'textbox', position, ...
               'String', label, 'EdgeColor', 'none', ...
               'Color', [0.10 0.10 0.10], ...
               'FontSize', font_size, 'FontWeight', 'normal', ...
               'Interpreter', 'tex');
end

function write_pdf(fig, out_path, resolution)
    out_dir = fileparts(out_path);
    if ~isempty(out_dir) && ~exist(out_dir, 'dir')
        mkdir(out_dir);
    end
    exportgraphics(fig, out_path, 'ContentType', 'image', 'Resolution', resolution);
    close(fig);
end

function require_file(path, hint)
    if ~exist(path, 'file')
        error('Missing %s. %s', path, hint);
    end
end

function label = time_label(seconds, unit)
    if ~isfinite(seconds)
        label = 't = n/a';
        return;
    end
    switch unit
        case 'as'
            label = sprintf('t = %.0f as', seconds / 1e-18);
        case 'fs'
            label = sprintf('t = %.0f fs', seconds / 1e-15);
        otherwise
            error('Unknown time label unit "%s".', unit);
    end
end

function label = trajectory_time_label(seconds)
    if seconds >= 1e-12
        label = sprintf('t = %.3g ps', seconds / 1e-12);
    else
        label = time_label(seconds, 'fs');
    end
end

function out_name = sphere_filename(seconds)
    as_value = round(seconds / 1e-18);
    out_name = sprintf('1s0_1fs-cloud_partial_%das.pdf', as_value);
end

function seconds = parse_time_values(t_final)
    if isnumeric(t_final)
        seconds = double(t_final(:)');
        for i = 1:numel(seconds)
            if seconds(i) > 1e-12
                seconds(i) = seconds(i) * 1e-18;
            end
        end
    else
        values = cellstr(string(t_final));
        seconds = zeros(1, numel(values));
        for i = 1:numel(values)
            seconds(i) = parse_time_string(values{i});
        end
    end
    if any(~isfinite(seconds)) || any(seconds <= 0)
        error('t_final must contain positive finite times.');
    end
end

function seconds = parse_time_string(value)
    tokens = regexp(strtrim(value), '^([0-9]*\.?[0-9]+(?:[eE][-+]?[0-9]+)?)\s*(as|fs|ps|s)$', ...
                    'tokens', 'once');
    if isempty(tokens)
        error('Could not parse t_final "%s". Use values such as "200as", "0.2fs", or "2e-16s".', value);
    end
    amount = str2double(tokens{1});
    switch lower(tokens{2})
        case 'as'
            seconds = amount * 1e-18;
        case 'fs'
            seconds = amount * 1e-15;
        case 'ps'
            seconds = amount * 1e-12;
        case 's'
            seconds = amount;
    end
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
