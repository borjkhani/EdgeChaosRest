function figS4_arclength()
%FIGS4_ARCLENGTH  Pseudo-arclength continuation of the quiet fixed point
%   (Supplementary Fig. S4; Experiment 3c). Mean E rate of the fixed point
%   versus gE_tot along the branch; solid: stable (no Jacobian eigenvalue
%   with positive real part, tau_N = 100 ms, thick blue), unstable: thin red. The branch
%   reaches a maximum of gE_tot (fold, circle) and turns back; triangle:
%   oscillatory instability preceding the fold.
S = fig_style(); fh = figure('Visible', 'off', 'Color', 'w');
rd = fullfile(fileparts(mfilename('fullpath')), '..', 'results', 'exp3c_arclength');
cases = [1 0; 1 0.5; 3 0]; ttl = {'connectivity 1, f_N = 0', 'connectivity 1, f_N = 0.5', 'connectivity 3, f_N = 0'};
for k = 1:3
    L = load(fullfile(rd, sprintf('seed%d_fN%02d.mat', cases(k,1), round(100*cases(k,2)))));
    o = L.out; st = o(:,6) == 0; ax = subplot(1, 3, k); hold(ax, 'on');
    plot(ax, o(:,1), o(:,3), '-', 'Color', [0.85 0.1 0.1], 'LineWidth', 1.3);
    y = o(:,3); y(~st) = NaN;
    plot(ax, o(:,1), y, '-', 'Color', S.col(1,:), 'LineWidth', 2.6);
    [~, i1] = max(o(1:find(diff(o(:,1)) < 0, 1), 1));
    h1 = plot(ax, o(i1,1), o(i1,3), 'o', 'Color', S.col(2,:), 'MarkerFaceColor', S.col(2,:), 'MarkerSize', S.ms + 1);
    iu = find(o(:,6) > 0, 1); h2 = [];
    if iu < i1, h2 = plot(ax, o(iu-1,1), o(iu-1,3), '^', 'Color', S.col(3,:), 'MarkerFaceColor', S.col(3,:), 'MarkerSize', S.ms + 1); end
    if k == 2, legend(ax, [h1 h2], {'first fold', 'oscillatory instability'}, 'Location', 'northwest', 'Box', 'off'); end
    xl = [1 max(o(o(:,3) < 6, 1)) + 0.1]; xlim(ax, xl); ylim(ax, [0 6]);
    xlabel(ax, 'g_E^{tot} (\times g_L)'); ylabel(ax, 'fixed-point E rate (Hz)');
    style_axes(ax, sprintf('%s  %s', char('a' + k - 1), ttl{k}));
end
save_fig(fh, 'FigS4_arclength', 14, 6.6);
close(fh);
end
