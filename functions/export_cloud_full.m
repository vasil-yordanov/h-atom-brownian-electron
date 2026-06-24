function export_cloud_full(mat_path, params_set_name, npts_cap, view_az, view_el, out_path)
    % EXPORT_CLOUD_FULL  Render a 3D position cloud (manuscript Fig. 1/2 style) from
    % the first NPTS_CAP trajectory points of a run -- every point, no decimation --
    % and write it to OUT_PATH as a rasterized PDF. This is the single renderer used
    % for all manuscript clouds (Figs. 1 and 2).
    %
    % Drawing the full 2p0 cloud (2e8 points) the naive way (double precision +
    % plot_trajectory_3d's copies) needs ~15-19 GB and OOM-kills on a 16 GB machine;
    % this reads each coordinate via matfile, casts to single, and frees
    % intermediates, peaking near ~5 GB. The bottleneck is then rasterizing the
    % segments (~12 min for the 2e8-point 2p0 cloud; seconds for the small 1s clouds).
    %
    %   mat_path        : run .mat holding r_traj/theta_traj/phi_traj, n, l (and optionally a_0)
    %   params_set_name : run label, for the cube size (cloud_cube_limit) and history range
    %   npts_cap        : draw the first min(npts_cap, N) points (Inf = the whole run)
    %   view_az,view_el : camera azimuth / elevation
    %   out_path        : output PDF path

    mf  = matfile(mat_path);
    n   = mf.n;
    l   = mf.l;
    a_0 = 5.29177e-11;
    try
        a_0 = double(mf.a_0);   % use the run's stored value if present
    catch
    end

    sz   = size(mf, 'r_traj');
    N    = sz(2);
    npts = min(npts_cap, N);

    % Read coordinates as single precision -- the whole array, or a leading
    % contiguous slice -- to keep peak memory low (avoids the double-precision copies
    % that OOM-kill a naive full draw).
    if npts >= N
        r  = single(mf.r_traj);
        th = single(mf.theta_traj);
        ph = single(mf.phi_traj);
    else
        r  = single(mf.r_traj(1, 1:npts));
        th = single(mf.theta_traj(1, 1:npts));
        ph = single(mf.phi_traj(1, 1:npts));
    end

    a0s = single(a_0);
    sth = sin(th);
    x = r .* sth .* cos(ph) ./ a0s;
    y = r .* sth .* sin(ph) ./ a0s;
    clear sth;
    z = r .* cos(th) ./ a0s;
    clear r th ph;

    fig = figure('Position', [100, 100, 620, 620], 'Color', 'w');
    plot3(x, y, z, '-', 'LineWidth', 1.0, 'Color', [0.1, 0.35, 0.8]);
    clear x y z;

    xlabel('X', 'FontSize', 14); ylabel('Y', 'FontSize', 14); zlabel('Z', 'FontSize', 14);
    axis equal; axis vis3d; camup([0 0 1]); set(gca, 'FontSize', 12);

    lim = cloud_cube_limit(n, l, a_0, params_set_name);
    xlim([-lim lim]); ylim([-lim lim]); zlim([-lim lim]);
    grid on;
    view(view_az, view_el);

    % Standardise for export: true cube box (camva refits the zoom AFTER the limits),
    % a margin for the tick numbers, and larger fonts (the 620 px figure is scaled
    % down ~3x in the manuscript, so 12-14 pt would come out tiny).
    set(gca, 'Position', [0.13, 0.13, 0.84, 0.83]);
    daspect([1 1 1]); pbaspect([1 1 1]); camva('auto');
    fs = 20; set(gca, 'FontSize', fs);
    set([get(gca, 'XLabel'), get(gca, 'YLabel'), get(gca, 'ZLabel')], 'FontSize', fs);

    % Remove the interactive axes toolbar (export/zoom/pan/rotate/home icons) so it
    % is not baked into the top-right of the rasterized export under -batch.
    set(gca, 'Toolbar', []);

    exportgraphics(fig, out_path, 'ContentType', 'image');
    close(fig);
    fprintf('%s written (%d points)\n', out_path, npts);
end
