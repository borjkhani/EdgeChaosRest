function fig3_regimes()
%FIG3_REGIMES  Excitatory gain x NMDA fraction (Experiment 1, 10 connectivities).
%   a  phase diagram: fraction of connectivities in the active state (mean E
%      rate >= 4 Hz) for every coupling and NMDA fraction; white markers and
%      bars: coupling at which P(active) = 0.5 from a logistic fit, with 95%
%      cluster-bootstrap confidence interval over connectivities
%   b  phase diagram of mean memory capacity (all runs)
%   c-e  mean E rate, autonomous exponent and memory capacity, summarised
%      separately for quiet runs (filled markers, solid lines) and active runs
%      (open markers, dashed lines): mean +- s.d. over the runs in that state
%   f  memory capacity versus conditional exponent for all distinct runs of
%      Experiments 1-3 (filled: quiet, open: active); black: mean +- s.e.m.
%      in bins of lambda_cond
S = fig_style(); D = load_scan('exp1_regime'); fN = [0 0.25 0.5 0.75];
gE = unique(D.pt(:,1))'; act = D.rateE >= 4;
fh = figure('Visible', 'off', 'Color', 'w');
P = nan(4, numel(gE)); M = P; N = P;
for k = 1:4
    for j = 1:numel(gE)
        s = abs(D.pt(:,2) - fN(k)) < 1e-9 & abs(D.pt(:,1) - gE(j)) < 1e-9;
        P(k,j) = mean(act(s)); M(k,j) = mean(D.MC(s)); N(k,j) = sum(s);
    end
end
dg = gE(2) - gE(1);
% a: phase diagram, probability of the active state
ax = subplot(2, 3, 1); hold(ax, 'on');
imagesc(ax, gE, 1:4, P); set(ax, 'CLim', [0 1], 'YDir', 'normal');
colormap(ax, seq_map([1 1 1], [0.10 0.27 0.55]));
for k = 1:4, for j = 1:numel(gE)
    tc = S.ink; if P(k,j) > 0.55, tc = [1 1 1]; end
    text(ax, gE(j), k, sprintf('%d', round(P(k,j) * N(k,j))), 'HorizontalAlignment', 'center', ...
         'FontSize', S.fs - 2, 'Color', tc);
end, end
for k = 1:4
    s = abs(D.pt(:,2) - fN(k)) < 1e-9;
    [g50, ci] = logistic_edge(D.pt(s,1), act(s), D.pt(s,3));
    line(ax, ci, (k + 0.32) * [1 1], 'Color', [0.85 0.1 0.1], 'LineWidth', 1.4);
    plot(ax, g50, k + 0.32, 'v', 'Color', [0.85 0.1 0.1], 'MarkerFaceColor', [0.85 0.1 0.1], 'MarkerSize', S.ms);
end
set(ax, 'YTick', 1:4, 'YTickLabel', {'0', '0.25', '0.5', '0.75'});
xlim(ax, [gE(1) - dg/2, gE(end) + dg/2]); ylim(ax, [0.5 4.5]);
xlabel(ax, 'g_E^{tot} (\times g_L)'); ylabel(ax, 'NMDA fraction f_N');
style_axes(ax, 'a  active connectivities (of 10)'); grid(ax, 'off');
colorbar(ax);
% b: phase diagram, memory capacity
ax = subplot(2, 3, 2); hold(ax, 'on');
imagesc(ax, gE, 1:4, M); set(ax, 'CLim', [0 4.5], 'YDir', 'normal');
colormap(ax, seq_map([1 1 1], [0.05 0.40 0.30]));
for k = 1:4, for j = 1:numel(gE)
    tc = S.ink; if M(k,j) > 2.6, tc = [1 1 1]; end
    text(ax, gE(j), k, sprintf('%.1f', M(k,j)), 'HorizontalAlignment', 'center', 'FontSize', S.fs - 3, 'Color', tc);
end, end
set(ax, 'YTick', 1:4, 'YTickLabel', {'0', '0.25', '0.5', '0.75'});
xlim(ax, [gE(1) - dg/2, gE(end) + dg/2]); ylim(ax, [0.5 4.5]);
xlabel(ax, 'g_E^{tot} (\times g_L)'); ylabel(ax, 'NMDA fraction f_N');
style_axes(ax, 'b  mean memory capacity'); grid(ax, 'off');
colorbar(ax);
% c-e: quiet and active states summarised separately
lab = {'c  mean E rate (Hz)', 'd  \lambda_{aut} (s^{-1})', 'e  memory capacity MC'};
yc = {'rateE', 'lam_aut', 'MC'};
for p = 1:3
    ax = subplot(2, 3, p + 2); hold(ax, 'on');
    for k = 1:4
        sel = abs(D.pt(:,2) - fN(k)) < 1e-9; dx = 0.02 * (k - 2.5);
        for st = 0:1
            m = nan(size(gE)); sd = m;
            for j = 1:numel(gE)
                s = sel & abs(D.pt(:,1) - gE(j)) < 1e-9 & act == st;
                if any(s), y = D.(yc{p})(s); m(j) = mean(y); sd(j) = std(y, 1); end
            end
            ok = find(isfinite(m)); ls = '-'; mf = S.col(k,:); if st, ls = '--'; mf = 'w'; end
            for j = ok, line(ax, [gE(j) gE(j)] + dx, m(j) + [-1 1] * sd(j), 'Color', 0.45 * S.col(k,:) + 0.55, 'LineWidth', 0.9); end
            hp = plot(ax, gE(ok) + dx, m(ok), [ls S.mk{k}], 'Color', S.col(k,:), 'MarkerFaceColor', mf, ...
                      'MarkerSize', S.ms - 0.5, 'LineWidth', 1.2);
            if st == 0, h(k) = hp; end %#ok<AGROW>
        end
    end
    xlim(ax, [1.1 3.1]); xlabel(ax, 'g_E^{tot} (\times g_L)'); style_axes(ax, lab{p});
    if p == 1, refline_h(ax, 4); ylim(ax, [0 26]);
        legend(ax, h, {'f_N = 0', '0.25', '0.5', '0.75'}, 'Location', 'north', 'Orientation', 'horizontal', 'Box', 'off'); end
    if p == 2, refline_h(ax, 0); ylim(ax, [-30 10]); end
    if p == 3, ylim(ax, [0 5]); end
end
% f: memory versus conditional exponent, Exp. 1-3 (distinct runs)
ax = subplot(2, 3, 6); hold(ax, 'on');
[lc, MC, rate] = pooled_runs();
q = rate < 4; cq = S.col(1,:); ca = S.col(2,:);
h1 = plot(ax, lc(q), MC(q), 'o', 'Color', cq, 'MarkerFaceColor', cq, 'MarkerSize', S.ms - 2.5);
h2 = plot(ax, lc(~q), MC(~q), 'o', 'Color', ca, 'MarkerSize', S.ms - 2);
ed = [-35 -20 -15 -10 -7.5 -5 -2.5 0 2.5 5 12]; xc = []; mm = []; se = [];
for b = 1:numel(ed) - 1
    s = lc >= ed(b) & lc < ed(b+1);
    if sum(s) >= 5, xc(end+1) = mean(lc(s)); mm(end+1) = mean(MC(s)); se(end+1) = std(MC(s)) / sqrt(sum(s)); end %#ok<AGROW>
end
for b = 1:numel(xc), line(ax, xc(b) * [1 1], mm(b) + [-1 1] * se(b), 'Color', S.ink, 'LineWidth', 1.2); end
h3 = plot(ax, xc, mm, '-s', 'Color', S.ink, 'MarkerFaceColor', S.ink, 'MarkerSize', S.ms - 1, 'LineWidth', 1.4);
xlim(ax, [-35 12]); ylim(ax, [0 5]); refline_v(ax, 0);
xlabel(ax, '\lambda_{cond} (s^{-1})'); ylabel(ax, 'MC');
style_axes(ax, 'f  MC vs \lambda_{cond} (Exp. 1-3)');
legend(ax, [h1 h2 h3], {'quiet', 'active', 'bin mean'}, 'Location', 'southwest', 'Box', 'off');
save_fig(fh, 'Fig3_regimes', 17, 11.5);
close(fh);
end

function [lc, MC, rate] = pooled_runs()
% distinct runs of Experiments 1-3 (reference conditions of Exp. 2 and 3
% duplicate Exp. 1 and are counted once)
E1 = load_scan('exp1_regime'); E2 = load_scan('exp2_kinetics'); E3 = load_scan('exp3_magnesium');
k2 = E2.pt(:,2) ~= 100; k3 = abs(E3.pt(:,2) - 1) > 1e-9;
f = @(n) [E1.(n); E2.(n)(k2); E3.(n)(k3)];
lc = f('lam_cond'); MC = f('MC'); rate = f('rateE');
end

function c = seq_map(c0, c1)
t = linspace(0, 1, 64)'; c = (1 - t) * c0 + t * c1;
end
