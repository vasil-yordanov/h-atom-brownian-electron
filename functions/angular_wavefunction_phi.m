function Y_phi = angular_wavefunction_phi(phi, m)
    % Azimuthal wavefunction (spherical harmonics) for hydrogen atom
    if m == 0
        Y_phi = ones(size(phi));
    elseif m > 0
        Y_phi = sqrt(2) * sin(abs(m) * phi);
    else
        Y_phi = sqrt(2) * cos(m * phi);
    end
end