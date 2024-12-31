function [r_min, r_max] = radial_histogram_range(n, l, a_0, params_set_name)
    % Function to provide radial histogram range
    % Inputs:
    %   n, l          : Quantum numbers (principal, angular, magnetic)
    %   a_0            : Scaling factor (e.g., Bohr radius)
    % Outputs:
    %   r_min         : the min radius
    %   r_max         : the max radius

    % Validate inputs
    if nargin < 3
        error('All three parameters (n, l, a_0) must be provided.');
    end

    % Initialize default values
     % Define initial values based on quantum numbers
    if n == 1 && l == 0 % 1s orbital
        if strcmp(params_set_name, 'fig_1')
            r_min = 0 * a_0;
            r_max = 1.5 * a_0;
        else
            r_min = 0.1* a_0;
            r_max = 5 * a_0;
        end
    elseif n == 2 && l == 0 % 2s orbital
        r_min = 0.1* a_0;
        r_max = 5 * (3 + sqrt(5)) * a_0;
    elseif n == 2 && l == 1  % 2p orbital
        r_min = 0.1* a_0;
        r_max = 5 * 4 * a_0;
    elseif n == 3 && l == 0 % 3s orbital
        r_min = 0.1* a_0;
        r_max = 5 * 0.8 * a_0;
    elseif n == 3 && l == 1 % 3p orbital
        r_min = 0.1* a_0;
        r_max = 5 * 4 * a_0;
    elseif n == 3 && l == 2 % 3d orbital
        r_min = 0.1* a_0;
        r_max = 5 * 9 * a_0;
    else
        error('Invalid combination of quantum numbers (n, l).');
    end
end