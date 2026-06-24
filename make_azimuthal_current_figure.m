% MAKE_AZIMUTHAL_CURRENT_FIGURE  Static two-panel manuscript figure for the
% phase-driven azimuthal circulation of the (2,1,m) states, m = +1, 0, -1
% (manuscript Section 9.1, Fig. 9).
%
%   (a) running time-average of the azimuthal drift b_phi, against the
%       analytical value <b_phi> = m (3*pi/32) hbar/(m_e a_0);
%   (b) unwrapped azimuthal coordinate phi(t), against the analytical
%       winding rate <dphi/dt> = m hbar/(8 m_e a_0^2).
%
% The coordinate distributions of the (2,1,1) state (manuscript Fig. 5) are
% exported by h_atom_post_processing, in the same style as Figs. 3-4.
%
% Requires the three production runs (1 ps each):
%   params_set_name='2p_m1_1ps';  h_atom
%   params_set_name='2p_m0_1ps';  h_atom
%   params_set_name='2p_mn1_1ps'; h_atom
% Then run this script from the repository root:
%   make_azimuthal_current_figure
%
% Output: figures/2p_pm1_1ps-azimuthal_current.pdf
% Copy the output next to the manuscript .tex file, which includes it under
% exactly this name.

addpath('functions');
[hbar, m_e, a_0, ~, ~, ~] = constants();

sets   = {'2p_m1_1ps', '2p_m0_1ps', '2p_mn1_1ps'};
labels = {'$m=+1$', '$m=0$', '$m=-1$'};
colors = {[0.85, 0.33, 0.10], [0.35, 0.35, 0.35], [0.00, 0.45, 0.74]};

% Analytical references for the (2,1,+-1) states (manuscript Eqs. (9.2)-(9.3)):
%   <b_phi>   = m * (3*pi/32) * hbar/(m_e a_0)   [<1/r> = 1/(4 a_0), <1/sin> = 3*pi/8]
%   <dphi/dt> = m * hbar/(8 m_e a_0^2)           [<1/r^2> = 1/(12 a_0^2), <1/sin^2> = 3/2]
b_phi_ref  = (3 * pi / 32) * hbar / (m_e * a_0);
phidot_ref = hbar / (8 * m_e * a_0^2);

S = cell(1, 3);
for k = 1:3
    f = fullfile('data', sets{k}, [sets{k} '.mat']);
    if ~exist(f, 'file')
        error('Missing %s -- run h_atom with params_set_name=''%s'' first.', f, sets{k});
    end
    S{k} = load(f, 'time_arr', 'b_phi_avg_traj', 'phi_traj', 'dt');
end

% Panel (a) shows the full running average of each run; panel (b) shows the
% per-step trajectory window common to all runs (at most 100 fs).
T_a = min(cellfun(@(s) s.time_arr(end), S));
T_b = min(1e-13, min(cellfun(@(s) numel(s.phi_traj) * s.dt, S)));

fig = figure('Units', 'centimeters', 'Position', [2, 2, 26, 8], 'Color', 'w');
% Two manually-placed panels, [left bottom width height] in normalized figure
% units, with a small explicit horizontal gap (tiledlayout presets were either
% too tight or too wide). Panel (a) ends at x=0.34; panel (b) starts at x=0.44,
% leaving room for its y-label plus a small gap.

% --- (a) running time-average of b_phi ----------------------------------
axes('Position', [0.065, 0.14, 0.275, 0.80]); hold on;
yline(0, 'k--', 'LineWidth', 0.8, 'HandleVisibility', 'off');
yline( b_phi_ref, 'r--', 'LineWidth', 1.2, 'HandleVisibility', 'off');
yline(-b_phi_ref, '--', 'Color', colors{3}, 'LineWidth', 1.2, 'HandleVisibility', 'off');
for k = 1:3
    plot(S{k}.time_arr, S{k}.b_phi_avg_traj, 'LineWidth', 1.6, ...
         'Color', colors{k}, 'DisplayName', labels{k});
end
hold off; grid on; box on;
xlim([0, T_a]);
ylim([-1e6, 1e6]);   % pin the tick scale (manual axes auto-picked a busier one)
xlabel('Time (s)', 'FontSize', 12);
ylabel('$\overline{b_\phi}$ (m/s)', 'Interpreter', 'latex', 'FontSize', 13);
lg_a = legend('Interpreter', 'latex', 'Location', 'east', 'Box', 'off', 'FontSize', 10);

% --- (b) unwrapped azimuth ------------------------------------------------
% Panel (b) shows the azimuth as the number of full revolutions phi/(2*pi):
% the analytical winding rate is phidot_ref/(2*pi) ~ 0.82 revolutions/fs.
axes('Position', [0.44, 0.14, 0.535, 0.80]); hold on;
tt = linspace(0, T_b, 50);
plot(tt,  phidot_ref * tt / (2*pi), 'r--', 'LineWidth', 1.2, 'HandleVisibility', 'off');
plot(tt, -phidot_ref * tt / (2*pi), '--', 'Color', colors{3}, 'LineWidth', 1.2, 'HandleVisibility', 'off');
for k = 1:3
    phi_u = unwrap(S{k}.phi_traj);
    phi_u = phi_u - phi_u(1);
    t_loc = (1:numel(phi_u)) * S{k}.dt;
    last = find(t_loc <= T_b, 1, 'last');
    idx = unique(round(linspace(1, last, 4000)));   % light vector PDF
    % The m=0 azimuth is a driftless random walk with large phi jumps during
    % pole visits (phi is a degenerate coordinate there); draw it thinner and
    % lighter so it reads as the null reference, not as a competing signal.
    if k == 2
        lw = 1.0; col = [0.62, 0.62, 0.62];
    else
        lw = 1.5; col = colors{k};
    end
    plot(t_loc(idx), phi_u(idx) / (2*pi), 'LineWidth', lw, 'Color', col, ...
         'DisplayName', labels{k});
end
hold off; grid on; box on;
xlim([0, T_b]);
xlabel('Time (s)', 'FontSize', 12);
ylabel('$\phi/2\pi$ (revolutions)', 'Interpreter', 'latex', 'FontSize', 13);
legend('Interpreter', 'latex', 'Location', 'northwest', 'Box', 'off', 'FontSize', 10);

% Both tiles now exist, so the tiled layout is final. Nudge panel (a)'s legend up
% into the gap between the y=0 line and the +b_phi line, so the y=0 reference line
% no longer runs through the "m=0" entry. (Must be done here, not right after the
% legend is created: the single-tile layout at that point is not the final geometry.)
drawnow;
lg_a.Position(2) = lg_a.Position(2) + 0.12;

if ~exist('figures', 'dir'); mkdir('figures'); end
out = fullfile('figures', '2p_pm1_1ps-azimuthal_current.pdf');
exportgraphics(fig, out, 'ContentType', 'vector');
fprintf('Figure written to %s\n', out);

% Console summary (numbers quoted in manuscript Section 9.1)
fprintf('\nAnalytic <b_phi>   = +/-%.4e m/s\n', b_phi_ref);
fprintf('Analytic <dphi/dt> = +/-%.4e rad/s (%.3f rad/fs)\n\n', phidot_ref, phidot_ref * 1e-15);
for k = 1:3
    bb = S{k}.b_phi_avg_traj(end);
    phi_u = unwrap(S{k}.phi_traj);
    t_loc = (1:numel(phi_u)) * S{k}.dt;
    last = find(t_loc <= T_b, 1, 'last');
    slope = (phi_u(last) - phi_u(1)) / t_loc(last);
    fprintf('%-11s  <b_phi>(%.0f fs) = %+1.4e m/s   winding slope(%.0f fs) = %+1.4e rad/s\n', ...
            sets{k}, S{k}.time_arr(end) * 1e15, bb, T_b * 1e15, slope);
end
