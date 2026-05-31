function ref = energy_reference_values(n, l, m)
    % ENERGY_REFERENCE_VALUES  Hard-coded analytical references for the
    % per-trajectory energy histograms. Fields are E_r, E_theta, E_phi, V, E
    % (all in joules). NaN means "no reference available for this state".

    ref.E_r = NaN; ref.E_theta = NaN; ref.E_phi = NaN;
    ref.V   = NaN; ref.E       = NaN;

    if n == 2 && l == 1 && abs(m) == 0       % (2,1,0): 2p_0
        ref.E_r     =  1.8166e-19;
        ref.E_theta =  3.6330e-19;            % angular KE total (m = 0 -> E_phi = 0)
        ref.E_phi   =  0;
        ref.V       = -10.899e-19;
        ref.E       = -5.4494e-19;
    elseif n == 2 && l == 0 && m == 0        % (2,0,0): 2s
        ref.E_r     =  5.4497e-19;
        ref.E_theta =  0;
        ref.E_phi   =  0;
        ref.V       = -10.899e-19;
        ref.E       = -5.4494e-19;
    elseif n == 1 && l == 0 && m == 0        % (1,0,0): 1s, virial theorem
        ref.E       = -2.1798e-18;
        ref.E_r     = -ref.E;
        ref.V       =  2 * ref.E;
        ref.E_theta =  0;
        ref.E_phi   =  0;
    end
end
