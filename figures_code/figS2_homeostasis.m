function figS2_homeostasis()
%FIGS2_HOMEOSTASIS  Robustness of homeostatic cycling to the scaling rule
%   (Supplementary Fig. S2; Experiment 8).
%   a-e  mean E rate (top of each panel band) and scaling factor theta
%        versus time after degeneration, connectivity 1, NMDA (fN = 0.5),
%        for five rules; dashed: homeostatic set point
%   f    fraction of time in the active state and time-averaged rate
%        relative to the set point, for all six networks and six rules
S = fig_style(); fh = figure('Visible', 'off', 'Color', 'w');
rd = fullfile(fileparts(mfilename('fullpath')), '..', 'results', 'exp8_homeostasis');
v = {'eta0.3', 'eta0.1', 'eta0.03', 'tanh', 'cont', 'floor'};
lab = {'\eta = 0.3 (reference)', '\eta = 0.1', '\eta = 0.03', 'bounded step (tanh, \eta = 0.1)', 'quasi-continuous (1-s epochs)', '\theta_{min} = 0.5'};
for k = 1:5
    L = load(fullfile(rd, sprintf('seed1_cond2_%s.mat', v{k})));
    t = (1:numel(L.H.rateE)) * L.h.T_epoch * 0.02 / 60;   % minutes
    ax = subplot(3, 2, k); hold(ax, 'on');
    h1 = plot(ax, t, L.H.rateE, '-', 'Color', S.col(2,:), 'LineWidth', 1);
    h2 = plot(ax, t, L.H.theta, '-', 'Color', S.col(1,:), 'LineWidth', 1.4);
    line(ax, [0 t(end)], L.target * [1 1], 'Color', S.ink2, 'LineStyle', '--', 'LineWidth', 0.8);
    set(ax, 'YScale', 'log'); ylim(ax, [0.2 30]); xlim(ax, [0 t(end)]);
    xlabel(ax, 'time after degeneration (min)'); style_axes(ax, sprintf('%s  %s', char('a' + k - 1), lab{k}));
    if k == 1, legend(ax, [h1 h2], {'mean E rate (Hz)', '\theta'}, 'Location', 'northeast', 'Box', 'off'); end
end
ax = subplot(3, 2, 6); hold(ax, 'on'); hh = [];
for k = 1:numel(v)
    for s = 1:3, for c = 1:2
        fn = fullfile(rd, sprintf('seed%d_cond%d_%s.mat', s, c, v{k}));
        if ~exist(fn, 'file'), continue; end
        L = load(fn); n = numel(L.H.rateE); i0 = round(n / 3) + 1;   % discard the first third
        fa = mean(L.H.rateE(i0:end) > 4); ra = mean(L.H.rateE(i0:end)) / L.target;
        mf = S.col(c,:); if c == 1, mf = 'w'; end
        plot(ax, k + 0.1 * (s - 2), fa, S.mk{c}, 'Color', S.col(c,:), 'MarkerFaceColor', mf, 'MarkerSize', S.ms);
    end, end
end
hh(1) = plot(ax, NaN, NaN, S.mk{1}, 'Color', S.col(1,:), 'MarkerSize', S.ms, 'LineWidth', 1.2);
hh(2) = plot(ax, NaN, NaN, S.mk{2}, 'Color', S.col(2,:), 'MarkerFaceColor', S.col(2,:), 'MarkerSize', S.ms);
set(ax, 'XTick', 1:numel(v), 'XTickLabel', {'0.3', '0.1', '0.03', 'tanh', 'cont', 'floor'});
xlim(ax, [0.5 numel(v) + 0.5]); ylim(ax, [0 1]); ylabel(ax, 'fraction of time active');
style_axes(ax, 'f  all networks (last 2/3 of each run)');
legend(ax, hh, {'AMPA only', 'NMDA'}, 'Location', 'northeast', 'Box', 'off');
save_fig(fh, 'FigS2_homeostasis', 14, 13);
close(fh);
end
