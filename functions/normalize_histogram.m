function hist_counts = normalize_histogram(hist_counts, hist_bins)
    bin_widths = diff(hist_bins);
    total_area = sum(hist_counts .* bin_widths);

    % Normalize the histogram counts so that the area under the histogram is 1
    if total_area > 0
        hist_counts = hist_counts / total_area;
    else
        hist_counts = hist_counts;  % Avoid division by zero
    end
end
