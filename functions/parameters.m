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

    scan_vmax_over_c = parse_scan_vmax_over_c(parameters_set);
    scan_dt = parse_scan_dt(parameters_set);
    if ~isnan(scan_vmax_over_c)
        M = 100;
        v_max_over_c = scan_vmax_over_c;
        if strcmp(parameters_set, 'test_scan_vmax_1p0c')
            n=2; l=1; m=0;
            n_steps = 1e7;
            dt = 1e-20;
            traj_points = 1;
        elseif ~isempty(regexp(parameters_set, '^2p0_scan_vmax_', 'once'))
            n=2; l=1; m=0;
            n_steps = 1e9;
            dt = 1e-20;
            traj_points = 1;
        elseif ~isempty(regexp(parameters_set, '^2s0_scan_vmax_', 'once'))
            n=2; l=0; m=0;
            n_steps = 1e9;
            dt = 1e-20;
            traj_points = 1;
        else
            error('The scan parameter set "%s" is not defined.', parameters_set);
        end
    elseif ~isnan(scan_dt)
        % Time-step (dt) convergence sweep. Like the vmax scan it uses
        % M=100 trajectories and traj_points=1 (statistical energy test, not
        % a trajectory-visualisation run), with v_max held at the production
        % cutoff. The total simulated time (n_steps*dt) is held fixed so all
        % dt points cover the same physical duration (equal physical time);
        % n_steps is rescaled accordingly. The dt value is encoded in the set
        % name, e.g. 2p0_scan_dt_1em20 -> dt = 1e-20, 5em21 -> 5e-21.
        M = 100;
        v_max_over_c = 1.0;
        traj_points = 1;
        dt = scan_dt;
        if ~isempty(regexp(parameters_set, '^test_scan_dt_', 'once'))
            n=2; l=1; m=0;
            total_time = 1e-13;
        elseif ~isempty(regexp(parameters_set, '^2p0_scan_dt_', 'once'))
            n=2; l=1; m=0;
            total_time = 1e-11;
        elseif ~isempty(regexp(parameters_set, '^2s0_scan_dt_', 'once'))
            n=2; l=0; m=0;
            total_time = 1e-11;
        else
            error('The dt scan parameter set "%s" is not defined.', parameters_set);
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

function v_max_over_c = parse_scan_vmax_over_c(parameters_set)
    tokens = regexp(parameters_set, '_scan_vmax_([0-9]+p[0-9]+c)$', 'tokens', 'once');
    if isempty(tokens)
        v_max_over_c = NaN;
        return;
    end

    token = strrep(tokens{1}, 'c', '');
    token = strrep(token, 'p', '.');
    v_max_over_c = str2double(token);
end

function dt = parse_scan_dt(parameters_set)
    % Decode the dt encoded in a *_scan_dt_* set name into a numeric value.
    % Encoding: mantissa 'p' = decimal point, 'e' separates the exponent,
    % 'm'/'p' after 'e' = minus/plus sign of the exponent.
    %   1em20   -> 1e-20      5em21  -> 5e-21
    %   1p5em20 -> 1.5e-20    1ep03  -> 1e+03
    tokens = regexp(parameters_set, '_scan_dt_([0-9p]+)e(m|p)([0-9]+)$', ...
                    'tokens', 'once');
    if isempty(tokens)
        dt = NaN;
        return;
    end

    mantissa = str2double(strrep(tokens{1}, 'p', '.'));
    exponent_sign = 1;
    if strcmp(tokens{2}, 'm')
        exponent_sign = -1;
    end
    exponent = str2double(tokens{3});
    dt = mantissa * 10^(exponent_sign * exponent);
end
