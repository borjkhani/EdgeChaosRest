function style_axes(ax, ttl, fs)
%STYLE_AXES  Recessive axes: light grid, no box, small ticks, panel title.
S = fig_style(); if nargin > 2 && ~isempty(fs), S.fs = fs; end
set(ax, 'FontSize', S.fs, 'Box', 'off', 'TickDir', 'out', ...
        'XColor', S.ink2, 'YColor', S.ink2, 'GridColor', S.grid);
grid(ax, 'on');
if nargin > 1 && ~isempty(ttl)
    title(ax, ttl, 'FontSize', S.fs + 1, 'FontWeight', 'bold', 'Color', S.ink);
end
end
