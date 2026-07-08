function [P_analytic_r, P_theta_analytic, P_phi_analytic] = manuscript_analytic_distributions(S)
    R_analytic = radial_wavefunction(S.r_values, S.n, S.l);
    P_analytic_r = S.r_values.^2 .* abs(R_analytic).^2;

    Y_theta_analytic = angular_wavefunction_theta(S.theta_values, S.l, S.m);
    P_theta_analytic = 2 * pi * sin(S.theta_values) .* abs(Y_theta_analytic).^2;

    P_phi_analytic = ones(size(S.phi_values)) / (2 * pi);
end
