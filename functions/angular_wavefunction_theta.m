function Y_theta = angular_wavefunction_theta(theta, l, m)
    % Polar wavefunction (spherical harmonics) for hydrogen atom
    if l == 0 && m == 0
        Y_theta = sqrt(1 / (4 * pi)) * ones(size(theta));
    elseif l == 1 && m == 0
        Y_theta = 0.5 * sqrt(3 / pi) * cos(theta);
    elseif l == 1 && m == 1
        Y_theta = - 0.5 * sqrt(3 / (2 * pi)) * sin(theta);
    elseif l == 2 && m == 0
        Y_theta = 0.25 * sqrt(5 / pi) * (3 * cos(theta).^2 - 1);
    else
        error('Angular wavefunction not implemented for l=%d, m=%d', l, m);
    end
end