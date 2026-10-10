function h = plot_states(ax, D, xcol, gcol, gvals, ycol, S, rate_thr)
%PLOT_STATES  Every run as a marker: filled = quiet state (mean E rate below
%   rate_thr Hz), open = active state. Lines join the means of the quiet
%   runs of each group (only where quiet runs exist).
if nargin < 8, rate_thr = 4; end
x = D.pt(:, xcol); g = D.pt(:, gcol); y = D.(ycol); quiet = D.rateE < rate_thr;
xs = unique(x)'; h = zeros(numel(gvals), 1);
for k = 1:numel(gvals)
    sel = abs(g - gvals(k)) < 1e-9;
    dx = 0.02 * (k - (numel(gvals) + 1) / 2);
    q = sel & quiet; a = sel & ~quiet;
    plot(ax, x(a) + dx, y(a), S.mk{k}, 'Color', S.col(k,:), 'MarkerSize', S.ms, 'LineWidth', 1);
    m = arrayfun(@(xx) mean(y(q & abs(x - xx) < 1e-9)), xs);
    ok = isfinite(m);
    plot(ax, xs(ok) + dx, m(ok), '-', 'Color', S.col(k,:), 'LineWidth', S.lw);
    h(k) = plot(ax, [NaN; x(q) + dx], [NaN; y(q)], S.mk{k}, 'Color', S.col(k,:), ...
                'MarkerFaceColor', S.col(k,:), 'MarkerSize', S.ms);
end
end
