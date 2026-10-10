function h = refline_v(ax, x, varargin)
%REFLINE_V  Vertical reference line across the current y-limits.
yl = get(ax, 'YLim'); h = line(ax, [x x], yl, 'Color', [0.45 0.45 0.45], ...
    'LineStyle', '--', 'LineWidth', 0.8, varargin{:});
end
