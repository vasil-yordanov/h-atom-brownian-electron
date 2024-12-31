function [n, l, m, n_steps, dt, traj_points, sigma_r_factor] = parameters(parameters_set)
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
    if strcmp(parameters_set, '1s0_1fs')
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
    elseif strcmp(parameters_set, '2p1_1ps')
        n=2; l=1; m=1;
        n_steps = 1e8;
        dt = 1e-20;
        traj_points = 1e7;
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
        error('The parameter set "%s" is not defined.', fig_name);
    end
end

