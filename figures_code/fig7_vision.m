function fig7_vision()
%FIG7_VISION  Retinal degeneration, homeostatic scaling and restoration (exp. 4).
%   a,b  scaling factor theta and mean E rate across homeostatic epochs after
%        loss of 75% of the feedforward drive (seed 1; AMPA-only and NMDA
%        networks). Dashed: healthy set point (b) and theta = 1 (a).
%   c,d  memory capacity and mean E rate at successive stages, for three
%        connectivity seeds (markers: circle 1, square 2, triangle 3) and both
%        synaptic conditions (colours):
%        H healthy; D degenerated, after scaling; R full restoration with
%        theta frozen; +Mg 2 mM Mg2+ at restoration; +blk 50% NMDA block at
%        restoration; S after re-scaling under restored input.
%        (+Mg and +blk act only through NMDA receptors and therefore leave the
%        AMPA-only network unchanged.)
%   e,f  restoration from actual states of the cycle (Experiment 11): mean E
%        rate and memory capacity after full restoration versus the scaling
%        factor theta at that moment (epochs 41-59), all connectivities
%        (markers) and conditions (colours: AMPA only, NMDA, NMDA + 2 mM Mg2+,
%        NMDA + 50% block).
S = fig_style(); S.fs = 10; fh = figure('Visible', 'off', 'Color', 'w');
here = fileparts(mfilename('fullpath'));
ld = @(s, c) getfield(load(fullfile(here, '..', 'results', 'exp6_vision', sprintf('seed%d_cond%d.mat', s, c))), 'Out');
condlab = {'AMPA only', 'NMDA (f_N = 0.5)'};
ax1 = subplot(2, 3, 1); hold(ax1, 'on'); ax2 = subplot(2, 3, 2); hold(ax2, 'on');
for c = 1:2
    O = ld(1, c); ep = 1:numel(O.hist_dep.theta);
    h1(c) = plot(ax1, ep, O.hist_dep.theta, ['-' S.mk{c}], 'Color', S.col(c,:), ...
        'MarkerFaceColor', S.col(c,:), 'MarkerSize', S.ms - 1, 'LineWidth', 1.2); %#ok<AGROW>
    plot(ax2, ep, O.hist_dep.rateE, ['-' S.mk{c}], 'Color', S.col(c,:), ...
        'MarkerFaceColor', S.col(c,:), 'MarkerSize', S.ms - 1, 'LineWidth', 1.2);
end
xlim(ax1, [0 61]); xlim(ax2, [0 61]); ylim(ax1, [0 3.9]); refline_h(ax1, 1); refline_h(ax2, ld(1, 2).healthy.rateE);
xlabel(ax1, 'homeostatic epoch'); ylabel(ax1, '\theta'); style_axes(ax1, 'a  scaling after degeneration', S.fs);
xlabel(ax2, 'homeostatic epoch'); ylabel(ax2, 'mean E rate (Hz)'); style_axes(ax2, 'b  rate after degeneration', S.fs);
legend(ax1, h1, condlab, 'Location', 'northwest', 'Box', 'off');
stages = {'H', 'D', 'R', '+Mg', '+blk', 'S'};
ax3 = subplot(2, 3, 3); hold(ax3, 'on'); ax4 = subplot(2, 3, 4); hold(ax4, 'on');
for c = 1:2
    for s = 1:3
        O = ld(s, c); R = O.rest(2);
        Ms = [O.healthy, O.deg_after, R.frozen, R.high_mg, R.nmda_block, R.rescaled];
        x = (1:6) + 0.12 * (c - 1.5) + 0.05 * (s - 2);
        mk = S.mk{s};
        plot(ax3, x, [Ms.MC], ['-' mk], 'Color', S.col(c,:), 'MarkerFaceColor', S.col(c,:), ...
             'MarkerSize', S.ms, 'LineWidth', 0.8);
        hh(c) = plot(ax4, x, [Ms.rateE], ['-' mk], 'Color', S.col(c,:), 'MarkerFaceColor', S.col(c,:), ...
             'MarkerSize', S.ms, 'LineWidth', 0.8); %#ok<AGROW>
    end
end
for a = [ax3 ax4]
    set(a, 'XTick', 1:6, 'XTickLabel', stages); xlim(a, [0.5 6.5]);
end
ylim(ax3, [0 5]); ylabel(ax3, 'memory capacity MC'); style_axes(ax3, 'c  memory for restored input', S.fs);
ylabel(ax4, 'mean E rate (Hz)'); style_axes(ax4, 'd  rate', S.fs); ylim(ax4, [0 12]);
for c = 1:2, hl(c) = plot(ax4, NaN, NaN, '-', 'Color', S.col(c,:), 'LineWidth', S.lw); end %#ok<AGROW>
legend(ax4, hl, condlab, 'Location', 'northwest', 'Box', 'off');
rd = fullfile(here, '..', 'results', 'exp11_restore');
ax5 = subplot(2, 3, 5); hold(ax5, 'on'); ax6 = subplot(2, 3, 6); hold(ax6, 'on'); hr = zeros(1, 4);
for s = 1:3, for c = 1:2
    L = load(fullfile(rd, sprintf('seed%d_cond%d.mat', s, c))); nm = fieldnames(L.R);
    for q = 1:numel(nm)
        ci = c + q - 1; X = L.R.(nm{q}); th = [X.theta] + 0.025 * (ci - 2.5);
        mf = S.col(ci,:);
        hr(ci) = plot(ax5, th, [X.rateE], S.mk{s}, 'Color', S.col(ci,:), 'MarkerFaceColor', mf, 'MarkerSize', S.ms - 1);
        plot(ax6, th, [X.MC], S.mk{s}, 'Color', S.col(ci,:), 'MarkerFaceColor', mf, 'MarkerSize', S.ms - 1);
    end
end, end
refline_h(ax5, 4); ylim(ax5, [0 20]); ylim(ax6, [0 5]);
for a = [ax5 ax6], xlim(a, [0 4.2]); xlabel(a, '\theta at restoration'); end
ylabel(ax5, 'mean E rate (Hz)'); ylabel(ax6, 'memory capacity MC');
style_axes(ax5, 'e  restoration from cycle states', S.fs); style_axes(ax6, 'f  memory after restoration', S.fs);
% colour key given in the caption (legend would cover data)
save_fig(fh, 'Fig7_vision', 15.5, 10.5);
close(fh);
end
