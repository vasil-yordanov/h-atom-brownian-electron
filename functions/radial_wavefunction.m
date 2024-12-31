function R = radial_wavefunction(r, n, l)
    % Numerical radial wavefunction for hydrogen atom
    % Bohr radius
    a_0 = 5.29177e-11;

    if n == 1 && l == 0
        % 1s state
        R = 2 * (1 / a_0)^(3/2) * exp(-r / a_0);
    elseif n == 2 && l == 0
        % 2s state
        R = 2 * (1 / (2 * a_0))^(3/2) .* (1 - r / (2 * a_0)) .* exp(-r / (2 * a_0));
    elseif n == 2 && l == 1
        % 2p state
        R = (1 / sqrt(3)) * (1 / (2 * a_0))^(3/2) .* (r / a_0) .* exp(-r / (2 * a_0));
    elseif n == 3 && l == 0
        % 3s state
        R = 2 * (1 / (3 * a_0))^(3/2) .* (1 - (2 * r) / (3 * a_0) + (2 * r.^2) / (27 * a_0^2)) .* exp(-r / (3 * a_0));
    elseif n == 3 && l == 1
        % 3p state
        R = (4 * sqrt(2) / 3) * (1 / (3 * a_0))^(3/2) .* (r / a_0) .* (1 - r / (6 * a_0)) .* exp(-r / (3 * a_0));
    elseif n == 3 && l == 2
        % 3d state
        R = (2 * sqrt(2) / (27 * sqrt(5))) * (1 / (3 * a_0))^(3/2) .* (r / a_0).^2 .* exp(-r / (3 * a_0));
    else
        % For higher n and l, use general formulas involving Laguerre polynomials
        error('Radial wave function not implemented for n = %d, l = %d', n, l);
    end
end

