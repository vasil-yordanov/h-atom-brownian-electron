function [deviation] = compute_distributions_deviation(hist_bins, hist_counts, values, P_analytic)
    % Normalize the Histogram Counts
    hist_counts_norm = normalize_histogram(hist_counts, hist_bins);

    % Calculate bin centers for R histogram
    bin_centers = hist_bins(1:end-1) + diff(hist_bins)/2;       
    P_analytic_bin = interp1(values, P_analytic, bin_centers, 'linear', 'extrap');

    %percentage_diff_bins = abs((hist_counts_norm - P_analytic_bin) ./ P_analytic_bin) * 100;
    %percentage_diff_bins(isnan(percentage_diff_bins)) = 0; % Set NaN (from 0/0) to 0
    %percentage_diff_bins(isinf(percentage_diff_bins)) = 0; % Set Inf (division by near-zero) to 0

    %absolute_diff = abs(hist_counts_norm - P_analytic_bin);

    % Compute weighted difference
    %deviation = sum(absolute_diff) / sum(hist_counts_norm) * 100;


    %total_weight = sum(hist_counts_norm);
    %deviation = sum(percentage_diff_bins .* P_analytic_bin) / total_weight; % weighted_diff

    % Compute Deviations for R, theta, and phi Distributions
    %deviation = sqrt(sum((hist_counts_norm - P_analytic_bin) .^2) * mean(diff(bin_centers)));

    % Calculate the overlap area
    overlap_area = sum(min(hist_counts_norm, P_analytic_bin)); % Intersection of histograms
    total_area = sum(max(hist_counts_norm, P_analytic_bin));   % Union of histograms
    
    % Compute overlap percentage
    deviation = 100 - (overlap_area / total_area) * 100; % overlap_percentage

end


