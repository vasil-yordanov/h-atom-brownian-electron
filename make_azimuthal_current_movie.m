% MAKE_AZIMUTHAL_CURRENT_MOVIE  "Comet" animation of the phase-driven azimuthal
% circulation of the (2,1,m) states, m = +1, 0, -1.
%
% The electron trajectory is drawn as a thin, fading azure trace over a short time
% window, viewed straight down the +z axis (the xy-plane). For m = +1 it sweeps
% counter-clockwise, for m = -1 clockwise, and for m = 0 it shows no net sense of
% rotation -- a direct visualisation of the azimuthal current
% b_phi = hbar*m/(m_e r sin theta) (manuscript Section 9.1, Fig. 9). Only the
% top-down view is produced, since that is where the rotation sense is legible.
%
% Reuses the stored trajectories of the three 1 ps production runs (no
% re-simulation needed):
%   params_set_name='2p_m1_1ps';  h_atom
%   params_set_name='2p_m0_1ps';  h_atom
%   params_set_name='2p_mn1_1ps'; h_atom
%
% It writes videos/2p_pm1_azimuthal_current_3panel.mp4.
% Convert to WebP with ffmpeg:
%   ffmpeg -y -i videos/2p_pm1_azimuthal_current_3panel.mp4 -vf "fps=50/3,scale=1500:-1:flags=lanczos" -c:v libwebp -q:v 85 -compression_level 6 -loop 0 -an movies/2p_pm1_azimuthal_current_3panel.webp

addpath('functions');
[~, ~, a_0] = constants();

out_mp4 = fullfile('videos', '2p_pm1_azimuthal_current_3panel.mp4');

% Single top-down (xy-plane) 3-panel figure (m = +1 | 0 | -1). The rotation
% sense (CCW vs CW) is only legible looking straight down +z, so no 3D view; and
% the side-by-side panels are the point, so no single-panel variant.
render_comet({'2p_m1_1ps','2p_m0_1ps','2p_mn1_1ps'}, {'$m=+1$','$m=0$','$m=-1$'}, ...
             out_mp4, a_0);

fprintf('\nMP4 written: %s\n', out_mp4);


function render_comet(sets, labels, out_mp4, a_0)
    % Render one MP4: top-down (xy-plane) comet for each set.
    N_win    = 1.0e6;    % stored points in the window (~4 revolutions for m=+1)
    n_frames = 150;      % animation frames
    L_tail   = 600000;   % comet tail length, stored points (longer trace)
    n_tail   = 1200;     % plotted points along the tail (decimated)
    lim      = 10;       % axis half-range in units of a_0

    np = numel(sets);
    X = cell(1, np);  Y = cell(1, np);
    for k = 1:np
        p = fullfile('data', sets{k}, [sets{k} '.mat']);
        if ~exist(p, 'file')
            error('Missing %s -- run h_atom with params_set_name=''%s'' first.', p, sets{k});
        end
        mf  = matfile(p);                        % partial read: only the window
        sz  = size(mf, 'r_traj');                % stored as a row vector (1 x n)
        nwk = min(N_win, max(sz));
        if sz(1) == 1
            r  = mf.r_traj(1, 1:nwk);  th = mf.theta_traj(1, 1:nwk);  ph = mf.phi_traj(1, 1:nwk);
        else
            r  = mf.r_traj(1:nwk, 1);  th = mf.theta_traj(1:nwk, 1);  ph = mf.phi_traj(1:nwk, 1);
        end
        r = r(:); th = th(:); ph = ph(:);
        X{k} = r .* sin(th) .* cos(ph) / a_0;
        Y{k} = r .* sin(th) .* sin(ph) / a_0;
    end

    out_dir = fileparts(out_mp4);
    if ~exist(out_dir, 'dir'); mkdir(out_dir); end

    % Panels sized for legibility (axis labels and titles readable in the animation).
    fig = figure('Color', 'w', 'Visible', 'off', ...
                 'Units', 'pixels', 'Position', [100, 100, 520 * np, 560]);

    vw = VideoWriter(out_mp4, 'MPEG-4');
    vw.FrameRate = 100 / 6;   % 150 frames over 9 seconds
    vw.Quality = 100;
    open(vw);

    % Fixed decimation grid: decimate the trajectory ONCE so every frame draws
    % the same sample points (only the visible range changes). Re-sampling the
    % tail per frame would shift the vertices and make the trace "swim"/jump.
    stride = max(1, round(L_tail / n_tail));
    samp   = 1:stride:N_win;
    heads  = round(linspace(round(N_win / n_frames), N_win, n_frames));

    for fI = 1:n_frames
        h = heads(fI);
        clf(fig);
        for k = 1:np
            if np > 1, ax = subplot(1, np, k); else, ax = axes('Parent', fig); end
            hold(ax, 'on');
            % Thin, fading azure trace. The tail is drawn as a series of short
            % line segments whose alpha ramps from transparent (old) to opaque
            % (new), giving a clean fading comet.
            azure = [0, 0.4470, 0.7410];
            lo    = max(1, h - L_tail);
            ti    = samp(samp >= lo & samp <= h);   % fixed-grid samples in the tail
            if numel(ti) >= 2
                nseg  = 40;
                edges = round(linspace(1, numel(ti), nseg + 1));
                for s = 1:nseg
                    seg  = edges(s):edges(s + 1);
                    aval = (s / nseg)^0.8;      % gentle fade: older tail stays visible longer
                    plot(ax, X{k}(ti(seg)), Y{k}(ti(seg)), '-', ...
                         'Color', [azure, aval], 'LineWidth', 0.5);
                end
            end
            axis(ax, 'equal');  xlim(ax, [-lim, lim]);  ylim(ax, [-lim, lim]);
            % Axis style matched to the other trajectory videos.
            grid(ax, 'on');  box(ax, 'on');  set(ax, 'FontSize', 11);
            xlabel(ax, 'X', 'FontSize', 13);  ylabel(ax, 'Y', 'FontSize', 13);
            title(ax, labels{k}, 'Interpreter', 'latex', 'FontSize', 15);
            hold(ax, 'off');
        end
        frame = getframe(fig);
        if ~isempty(vw.Height) && ...
                (size(frame.cdata, 1) ~= vw.Height || size(frame.cdata, 2) ~= vw.Width)
            frame.cdata = imresize(frame.cdata, [vw.Height, vw.Width]);
        end
        writeVideo(vw, frame);
    end
    close(vw);
    close(fig);
    fprintf('%d frames -> %s\n', n_frames, out_mp4);
end
