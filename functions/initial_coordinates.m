function [r_init, theta_init, phi_init] = initial_coordinates(n, l, m, M, a_0, params_set_name)
    % Function to provide initial coordinates of the particles
    % Inputs:
    %   n, l, m       : Quantum numbers (principal, angular, magnetic)
    %   M             : Number of particles
    %   a_0            : Scaling factor (e.g., Bohr radius)
    % Outputs:
    %   r_init        : Initial radius values
    %   theta_init    : Initial polar angle values
    %   phi_init      : Initial azimuthal angle values

    % Validate inputs
    if nargin < 5
        error('All five parameters (n, l, m, M, a_0) must be provided.');
    end

    % Initialize default values
     % Define initial values based on quantum numbers
    if n == 1 && l == 0 && m == 0 % 1s0 orbital
        if strcmp(params_set_name, '1s0_1ps')
            r_init = 10 * a_0 * ones(M, 1);
        else
            r_init = a_0 * ones(M, 1);
        end
        theta_init = pi / 2 * ones(M, 1);
        phi_init = zeros(M, 1);
    elseif n == 2 && l == 0 && m == 0 % 2s0 orbital
        r_init = (3 + sqrt(5)) * a_0 * ones(M, 1);
        theta_init = pi / 2 * ones(M, 1);
        phi_init = zeros(M, 1);
    elseif n == 2 && l == 1 && m == 0 % 2p0 orbital
        r_init = 4 * a_0 * ones(M, 1);
        theta_init = pi / 3 * ones(M, 1);
        phi_init = zeros(M, 1);
    elseif n == 2 && l == 1 && m == 1 % 2p1 orbital
        r_init = 4 * a_0 * ones(M, 1);
        theta_init = pi / 3 * ones(M, 1);
        phi_init = zeros(M, 1);
    elseif n == 3 && l == 0 && m == 0 % 3s0 orbital
        r_init = 0.8 * a_0 * ones(M, 1);
        theta_init = pi / 3 * ones(M, 1);
        phi_init = zeros(M, 1);
    elseif n == 3 && l == 1 && m == 0 % 3p0 orbital
        r_init = 4 * a_0 * ones(M, 1);
        theta_init = pi / 3 * ones(M, 1);
        phi_init = zeros(M, 1);
    elseif n == 3 && l == 2 && m == 0 % 3d0 orbital
        r_init = 9 * a_0 * ones(M, 1);
        theta_init = pi / 3 * ones(M, 1);
        phi_init = zeros(M, 1);
    else
        error('Invalid combination of quantum numbers (n, l, m).');
    end
end