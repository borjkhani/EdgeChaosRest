function h = refline_h(ax, y, varargin)
%REFLINE_H  Horizontal reference line across the current x-limits.
xl = get(ax, 'XLim'); h = line(ax, xl, [y y], 'Color', [0.45 0.45 0.45], ...
    'LineStyle', '--', 'LineWidth', 0.8, varargin{:});
end
