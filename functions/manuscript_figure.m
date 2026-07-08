function fig = manuscript_figure(default_font_size)
    if nargin < 1 || isempty(default_font_size)
        default_font_size = 20;
    end
    fig = figure('Visible', 'off', 'Color', 'w');
    set(fig, 'DefaultAxesFontSize', default_font_size, ...
             'DefaultTextFontSize', default_font_size, ...
             'DefaultLegendFontSize', default_font_size);
end
