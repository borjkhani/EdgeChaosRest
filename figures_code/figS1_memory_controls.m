function figS1_memory_controls()
%FIGS1_MEMORY_CONTROLS  Controls for the memory results (Supplementary Fig. S1).
%   a  memory capacity versus distance to each network's own linear-stability
%      edge, d = 1 - gE/g_edge (continuation, Exp. 3b), for connectivities
%      1-3 of Exp. 1 (gain x fN), Exp. 2 (gain x tau_N) and Exp. 3 (gain x Mg)
%   b  delay-resolved memory MC_k, mean +- s.e.m. over quiet and active runs
%      of Exp. 1 (10 connectivities)
%   c  quadratic (nonlinear) memory capacity versus lambda_cond (Exp. 1);
%      dashed: chance level for an input never presented
%   d  linear and quadratic capacity: quiet and active runs of Exp. 1,
%      networks without recurrence (same drive; drive matched to the quiet
%      rate) and chance level; mean +- s.d.
S = fig_style(); fh = figure('Visible', 'off', 'Color', 'w');
[d, mc, q, grp] = distance_data();
ax = subplot(2, 2, 1); hold(ax, 'on'); hs = [];
for k = 1:3
    s = grp == k;
    plot(ax, d(s & ~q), mc(s & ~q), S.mk{k}, 'Color', S.col(k,:), 'MarkerSize', S.ms - 1);
    hs(k) = plot(ax, [NaN; d(s & q)], [NaN; mc(s & q)], S.mk{k}, 'Color', S.col(k,:), 'MarkerFaceColor', S.col(k,:), 'MarkerSize', S.ms - 1); %#ok<AGROW>
end
refline_v(ax, 0); xlim(ax, [-1.2 0.5]); ylim(ax, [0 5]);
xlabel(ax, 'distance to own edge, 1 - g_E^{tot}/g_{edge}'); ylabel(ax, 'MC');
style_axes(ax, 'a  MC vs distance to the linear edge');
legend(ax, hs, {'gain x f_N', 'gain x \tau_N', 'gain x [Mg^{2+}]'}, 'Location', 'northwest', 'Box', 'off');
% b
D = load_scan('exp1_regime'); qa = D.rateE < 4; K = size(D.MCk, 2); lag = (1:K) * 20;
ax = subplot(2, 2, 2); hold(ax, 'on'); cc = {S.col(1,:), S.col(2,:)}; st = {qa, ~qa}; hb = [];
for k = 1:2
    m = mean(D.MCk(st{k}, :), 1); se = std(D.MCk(st{k}, :), 0, 1) / sqrt(sum(st{k}));
    for j = 1:K, line(ax, lag(j) * [1 1], m(j) + [-1 1] * se(j), 'Color', cc{k}, 'LineWidth', 0.8); end
    hb(k) = plot(ax, lag, m, ['-' S.mk{k}], 'Color', cc{k}, 'MarkerFaceColor', cc{k}, 'MarkerSize', S.ms - 2, 'LineWidth', S.lw); %#ok<AGROW>
end
xlim(ax, [0 400]); ylim(ax, [0 1]); xlabel(ax, 'lag (ms)'); ylabel(ax, 'MC_k');
style_axes(ax, 'b  delay-resolved memory (Exp. 1)');
legend(ax, hb, {sprintf('quiet (n=%d)', sum(qa)), sprintf('active (n=%d)', sum(~qa))}, 'Location', 'northeast', 'Box', 'off');
% c
E = load_scan_extra('exp10_nonlinear');
ax = subplot(2, 2, 3); hold(ax, 'on'); qe = E.rateE < 4;
h1 = plot(ax, E.lam_cond(qe), E.NLC(qe), 'o', 'Color', S.col(1,:), 'MarkerFaceColor', S.col(1,:), 'MarkerSize', S.ms - 2.5);
h2 = plot(ax, E.lam_cond(~qe), E.NLC(~qe), 'o', 'Color', S.col(2,:), 'MarkerSize', S.ms - 2);
line(ax, [-35 12], mean(E.NLC_shuf) * [1 1], 'Color', S.ink2, 'LineStyle', '--', 'LineWidth', 0.9);
refline_v(ax, 0); xlim(ax, [-35 12]); xlabel(ax, '\lambda_{cond} (s^{-1})'); ylabel(ax, 'quadratic capacity');
style_axes(ax, 'c  nonlinear memory vs \lambda_{cond}');
legend(ax, [h1 h2], {'quiet', 'active'}, 'Location', 'northeast', 'Box', 'off');
% d
L = load(fullfile(fileparts(mfilename('fullpath')), '..', 'results', 'exp9_numerics', 'numerics.mat'));
sd0 = arrayfun(@(b) b.same_drive.MC, L.base); sm0 = arrayfun(@(b) b.matched.MC, L.base);
sd1 = arrayfun(@(b) b.same_drive.NLC, L.base); sm1 = arrayfun(@(b) b.matched.NLC, L.base);
vals = {E.MC(qe), E.MC(~qe), sd0, sm0, E.MC_shuf; E.NLC(qe), E.NLC(~qe), sd1, sm1, E.NLC_shuf};
ax = subplot(2, 2, 4); hold(ax, 'on'); hd = [];
for r = 1:2
    for c = 1:5
        v = vals{r, c}; x = c + 0.15 * (2 * r - 3);
        line(ax, x * [1 1], mean(v) + [-1 1] * std(v), 'Color', S.col(r,:), 'LineWidth', 1);
        hh = plot(ax, x, mean(v), S.mk{r}, 'Color', S.col(r,:), 'MarkerFaceColor', S.col(r,:), 'MarkerSize', S.ms);
        if c == 1, hd(r) = hh; end %#ok<AGROW>
    end
end
set(ax, 'XTick', 1:5, 'XTickLabel', {'quiet', 'active', 'no rec.', 'no rec.*', 'chance'});
xlim(ax, [0.4 5.6]); ylim(ax, [0 6]); ylabel(ax, 'capacity');
style_axes(ax, 'd  baselines');
legend(ax, hd, {'linear MC', 'quadratic'}, 'Location', 'northeast', 'Box', 'off');
save_fig(fh, 'FigS1_memory_controls', 14, 11.5);
close(fh);
end

function [d, mc, q, grp] = distance_data()
root = fullfile(fileparts(mfilename('fullpath')), '..');
rd = fullfile(root, 'results', 'exp3b_continuation');
ld = @(s, f, m) load(fullfile(rd, sprintf('seed%d_fN%02d_Mg%03d.mat', s, round(100*f), round(100*m))));
d = []; mc = []; q = []; grp = [];
E1 = load_scan('exp1_regime');
for s = 1:3, for f = [0 0.25 0.5 0.75]
    e = edge_of(ld(s, f, 1), 3); i = E1.pt(:,3) == s & E1.pt(:,2) == f;
    d = [d; 1 - E1.pt(i,1) / e]; mc = [mc; E1.MC(i)]; q = [q; E1.rateE(i) < 4]; grp = [grp; ones(sum(i), 1)]; %#ok<AGROW>
end, end
E2 = load_scan('exp2_kinetics'); tN = [25 50 100 200];
for s = 1:2, for t = 1:4
    if tN(t) == 100, continue; end
    e = edge_of(ld(s, 0.5, 1), t); i = E2.pt(:,3) == s & E2.pt(:,2) == tN(t);
    d = [d; 1 - E2.pt(i,1) / e]; mc = [mc; E2.MC(i)]; q = [q; E2.rateE(i) < 4]; grp = [grp; 2 * ones(sum(i), 1)]; %#ok<AGROW>
end, end
E3 = load_scan('exp3_magnesium');
for s = 1:2, for m = [0.5 2]
    e = edge_of(ld(s, 0.5, m), 3); i = E3.pt(:,3) == s & abs(E3.pt(:,2) - m) < 1e-9;
    d = [d; 1 - E3.pt(i,1) / e]; mc = [mc; E3.MC(i)]; q = [q; E3.rateE(i) < 4]; grp = [grp; 3 * ones(sum(i), 1)]; %#ok<AGROW>
end, end
q = logical(q);
end

function e = edge_of(L, t)
if isfinite(L.gE_unstable(t)), e = L.gE_unstable(t);
elseif isfield(L, 'gE_end') && isfinite(L.gE_end), e = L.gE_end;
else, e = L.gEs(end); end
end

function D = load_scan_extra(name)
here = fileparts(mfilename('fullpath')); old = cd(fullfile(here, '..'));
T = collect_grid(name); cd(old);
D.pt = reshape([T.pt], 3, [])';
for f = {'MC', 'MC_shuf', 'NLC', 'NLC_shuf', 'lam_cond', 'rateE'}
    D.(f{1}) = arrayfun(@(t) double(t.(f{1})), T)';
end
end
