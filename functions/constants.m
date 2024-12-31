function [hbar ,m_e, a_0, e_charge, epsilon_0, c] = constants()
    %CONSTANTS Summary of this function goes here
    %   Detailed explanation goes here

    hbar = 1.0545718e-34;        % Reduced Planck constant (J*s)
    m_e = 9.10938356e-31;        % Mass of an electron (kg)
    a_0 = 5.29177e-11;           % Bohr radius (m)
    e_charge = 1.602176634e-19;  % Elementary charge (C)
    epsilon_0 = 8.854187817e-12; % Vacuum permittivity (F/m)
    c=3e8;                       % speed of light (m/s)
end

