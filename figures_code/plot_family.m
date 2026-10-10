function h = plot_family(ax, D, xcol, gcol, gvals, ycol, S, transform)
%PLOT_FAMILY  Seed-mean line per group value plus individual seeds as dots.
%   x = D.pt(:,xcol), group = D.pt(:,gcol); y = D.(ycol).
if nargin < 8, transform = @(y) y; end
h = zeros(numel(gvals), 1);
x = D.pt(:, xcol); g = D.pt(:, gcol); y = transform(D.(ycol));
xs = unique(x)';
for k = 1:numel(gvals)
    sel = abs(g - gvals(k)) < 1e-9;
    m = arrayfun(@(xx) mean(y(sel & abs(x - xx) < 1e-9), 'omitnan'), xs);
    plot(ax, x(sel) + 0.012 * (k - 2.5), y(sel), '.', 'Color', 0.55 * S.col(k,:) + 0.45, 'MarkerSize', 7);
    h(k) = plot(ax, xs, m, ['-' S.mk{k}], 'Color', S.col(k,:), 'LineWidth', S.lw, ...
                'MarkerFaceColor', S.col(k,:), 'MarkerSize', S.ms);
end
end
