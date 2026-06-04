function [n, l, m, n_steps, dt, traj_points, sigma_r_factor, M, v_max_over_c] = parameters(parameters_set)
    % PARAMETERS returns some main parameters for the simulation
    % depending on the input name.
    %
    % Usage:
    %   [n, l, m, n_steps, dt] = parameters('fig_1');
    %
    % Input parameter parameters_set can be: 
    %   - 'fig_1'
    %   - 'fig_2a'
    %   - 'fig_2b'
    %   - 'fig_3'
    %   - 'fig_4_6b_7b'
    %   - 'fig_5_6a_7a'
    %
    % Output parameters are:
    % n, l, m     - Quantum numbers
    % n_steps     - Total number of simulation steps
    % dt          - Time step.
    % traj_points - Number of particle coordinates to store

    sigma_r_factor = 1;
    M = 1;
    v_max_over_c = 1.0;

    [scan_vmax_over_c, scan_dt] = parse_scan(parameters_set);
    if ~isnan(scan_vmax_over_c)
        % Unified kinetic scan. The set name encodes both v_max/c and the
        % integration step dt in zeptoseconds (1 zs = 1e-21 s); the dt tag is
        % required (e.g. 2s0_scan_vmax_3p0c_dt_1zs). Both v_max/c and dt can be
        % overridden at run time (vmax_override / dt_override). M = 100
        % trajectories, traj_points = 1 (a statistical energy test). The total
        % simulated time is fixed per state and n_steps = round(total_time/dt),
        % so every dt covers the same physical duration (equal physical time).
        M = 100;
        v_max_over_c = scan_vmax_over_c;
        dt = scan_dt;
        traj_points = 1;
        if ~isempty(regexp(parameters_set, '^test_scan_vmax_', 'once'))
            n=2; l=1; m=0;
            total_time = 1e-13;
        elseif ~isempty(regexp(parameters_set, '^2p0_scan_vmax_', 'once'))
            n=2; l=1; m=0;
            total_time = 1e-11;
        elseif ~isempty(regexp(parameters_set, '^2s0_scan_vmax_', 'once'))
            n=2; l=0; m=0;
            total_time = 1e-11;
        else
            error('The scan parameter set "%s" is not defined.', parameters_set);
        end
        n_steps = round(total_time / dt);
    elseif strcmp(parameters_set, '1s0_1fs')
        n=1; l=0; m=0;
        n_steps = 1e5;
        dt = 1e-20;
        traj_points = 1e5;
        sigma_r_factor = 0;
    elseif strcmp(parameters_set, '1s0_1ps') 
        n=1; l=0; m=0;
        n_steps = 1e8;
        dt = 1e-20;
        traj_points = 1e7;
    elseif strcmp(parameters_set, '2p0_10ps')
        n=2; l=1; m=0;
        n_steps = 1e9;
        dt = 1e-20;
        traj_points = 2e8;
    elseif strcmp(parameters_set, '2p_m0_1ps') || strcmp(parameters_set, '2p0_1ps')
        n=2; l=1; m=0;
        n_steps = 5e6;
        dt = 10e-21;
        traj_points = 2e6;
    elseif strcmp(parameters_set, '2p_m1_1ps') || strcmp(parameters_set, '2p1_1ps')
        n=2; l=1; m=1;
        n_steps = 2e6;
        dt = 5e-21;
        traj_points = 2e6;
    elseif strcmp(parameters_set, '2p_mn1_1ps')
        n=2; l=1; m=-1;
        n_steps = 2e6;
        dt = 5e-21;
        traj_points = 2e6;
    elseif strcmp(parameters_set, '2s0_10ps')
        n=2; l=0; m=0;
        n_steps = 1e9;
        dt = 1e-20;
        traj_points = 1e7;
    elseif strcmp(parameters_set, 'test')
        n=1; l=0; m=0;
        n_steps = 1e5;
        dt = 1e-20;
        traj_points = 1e5;
        sigma_r_factor = 0;
    else
        error('The parameter set "%s" is not defined.', parameters_set);
    end
end

function [v_max_over_c, dt] = parse_scan(parameters_set)
    % Parse v_max/c and dt from a unified kinetic-scan set name, e.g.
    %   2s0_scan_vmax_3p0c_dt_1zs    -> v_max/c = 3.0, dt = 1e-21 s
    %   2s0_scan_vmax_3p0c_dt_10zs   -> v_max/c = 3.0, dt = 1e-20 s
    % v_max/c uses 'p' for the decimal point and a trailing 'c'. dt is encoded
    % in zeptoseconds (1 zs = 1e-21 s), again with 'p' for the decimal point,
    % and is REQUIRED: a scan name without a _dt_<d>zs tag is an error.
    % Returns v_max/c = NaN (and dt = NaN) when the name is not a scan set.
    v_tokens = regexp(parameters_set, '_scan_vmax_([0-9]+p[0-9]+c)', ...
                      'tokens', 'once');
    if isempty(v_tokens)
        v_max_over_c = NaN;
        dt = NaN;
        return;
    end
    token = strrep(strrep(v_tokens{1}, 'c', ''), 'p', '.');
    v_max_over_c = str2double(token);

    dt_tokens = regexp(parameters_set, '_dt_([0-9p]+)zs', 'tokens', 'once');
    if isempty(dt_tokens)
        error(['Scan set name "%s" must include an explicit dt tag, e.g. ' ...
               '"%s_dt_1zs" (dt in zeptoseconds, 1 zs = 1e-21 s).'], ...
              parameters_set, parameters_set);
    end
    dt = str2double(strrep(dt_tokens{1}, 'p', '.')) * 1e-21;
end
