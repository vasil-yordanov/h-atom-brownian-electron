function manuscript_bump_fonts(fs)
    if nargin < 1 || isempty(fs)
        fs = 20;
    end
    set(gca, 'FontSize', fs);
    set([get(gca, 'XLabel'), get(gca, 'YLabel')], 'FontSize', fs + 2);
    lg = get(gca, 'Legend');
    if ~isempty(lg)
        lg.FontSize = fs - 2;
    end
end
