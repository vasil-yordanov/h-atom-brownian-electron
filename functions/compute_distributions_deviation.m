function [deviation] = compute_distributions_deviation(hist_bins, hist_counts, values, P_analytic)
    % Normalize the Histogram Counts
    hist_counts_norm = normalize_histogram(hist_counts, hist_bins);

    % Calculate bin centers for R histogram
    bin_centers = hist_bins(1:end-1) + diff(hist_bins)/2;       
    P_analytic_bin = interp1(values, P_analytic, bin_centers, 'linear', 'extrap');

    % Calculate the overlap area
    overlap_area = sum(min(hist_counts_norm, P_analytic_bin)); % Intersection of histograms
    total_area = sum(max(hist_counts_norm, P_analytic_bin));   % Union of histograms
    
    % Compute overlap percentage
    deviation = 100 - (overlap_area / total_area) * 100; % overlap_percentage

end


