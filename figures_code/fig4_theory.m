function fig4_theory()
%FIG4_THEORY  How the quiet state loses stability (Newton continuation).
%   a  static gain rho_K along the quiet fixed-point branch (connectivity 1,
%      four NMDA fractions). rho_K and the branch do not depend on tau_N.
%      Open markers: oscillatory (Hopf) instability for tau_N = 100 ms; the
%      branch ends at a fold.
%   b  largest real part of the Jacobian eigenvalues along the branch for
%      four NMDA decay constants (connectivity 1, fN = 0.5); the zero
%      crossings (Hopf points) move with tau_N.
%   c  coupling at which the quiet state is lost (first instability: Hopf if
%      it occurs before the fold, otherwise the fold) versus tau_N, for three
%      connectivities (markers) and three NMDA fractions (colours); dashed:
%      fold of the AMPA-only network of the same connectivity.
%   d  the same versus [Mg2+] (fN = 0.5, tau_N = 100 ms).
%   e  autonomous exponent of quiet fixed points in Experiment 2; dashed
%      lines: -1/tau_N.
%   f  rho_K with versus without the NMDA negative slope conductance
%      (Experiment 3).
S = fig_style(); fh = figure('Visible', 'off', 'Color', 'w');
% larger text, lines and markers than the other figures (dense panels)
S.fs = 10; S.lw = 2; S.ms = 6; sa = @(ax, t) style_axes(ax, t, S.fs);
here = fileparts(mfilename('fullpath')); rd = fullfile(here, '..', 'results', 'exp3b_continuation');
ld = @(s, f, m) load(fullfile(rd, sprintf('seed%d_fN%02d_Mg%03d.mat', s, round(100*f), round(100*m))));
edge = @(L, t) pick_edge(L, t);
fNs = [0 0.25 0.5 0.75]; tauN = [25 50 100 200];
% a
ax = subplot(2, 3, 1); hold(ax, 'on');
for k = 1:4
    L = ld(1, fNs(k), 1);
    h(k) = plot(ax, L.gEs, L.rhoK, '-', 'Color', S.col(k,:), 'LineWidth', S.lw); %#ok<AGROW>
    plot(ax, L.gEs(end), L.rhoK(end), S.mk{k}, 'Color', S.col(k,:), 'MarkerFaceColor', S.col(k,:), 'MarkerSize', S.ms);
    if isfinite(L.gE_unstable(3))
        j = find(abs(L.gEs - L.gE_unstable(3)) < 1e-9, 1);
        plot(ax, L.gEs(j), L.rhoK(j), S.mk{k}, 'Color', S.col(k,:), 'MarkerSize', S.ms + 2, 'LineWidth', 1.2);
    end
end
xlim(ax, [1 2.3]); ylim(ax, [0.4 1.4]); refline_h(ax, 1);
xlabel(ax, 'g_E^{tot} (\times g_L)', 'FontSize', S.fs); ylabel(ax, '\rho_K', 'FontSize', S.fs);
sa(ax, 'a  static gain on the quiet branch');
legend(ax, h, {'f_N = 0', '0.25', '0.5', '0.75'}, 'Location', 'southeast', 'Box', 'off');
% b
ax = subplot(2, 3, 2); hold(ax, 'on'); L = ld(1, 0.5, 1);
for t = 1:4
    hb(t) = plot(ax, L.gEs, L.maxRe(:, t), '-', 'Color', S.col(t,:), 'LineWidth', S.lw); %#ok<AGROW>
    if isfinite(L.gE_unstable(t)), plot(ax, L.gE_unstable(t), 0, S.mk{t}, 'Color', S.col(t,:), 'MarkerFaceColor', S.col(t,:), 'MarkerSize', S.ms); end
end
xlim(ax, [1 2.2]); ylim(ax, [-25 8]); refline_h(ax, 0);
xlabel(ax, 'g_E^{tot} (\times g_L)', 'FontSize', S.fs); ylabel(ax, 'max Re \mu (s^{-1})', 'FontSize', S.fs);
sa(ax, 'b  oscillatory instability moves with \tau_N');
legend(ax, hb, {'\tau_N = 25 ms', '50', '100', '200'}, 'Location', 'northwest', 'Box', 'off');
% c
ax = subplot(2, 3, 3); hold(ax, 'on');
for s = 1:3
    L0 = ld(s, 0, 1);
    line(ax, [20 230], L0.gEs(end) * [1 1], 'Color', [0.35 0.35 0.35], 'LineStyle', '--', 'LineWidth', 1.1);
    for k = 2:4
        L = ld(s, fNs(k), 1); e = arrayfun(@(t) edge(L, t), 1:4);
        hp = plot(ax, tauN, e, ['-' S.mk{s}], 'Color', S.col(k,:), 'MarkerFaceColor', S.col(k,:), ...
             'MarkerSize', S.ms, 'LineWidth', 1.6);
        if s == 1, hc(k-1) = hp; end %#ok<AGROW>
    end
end
set(ax, 'XScale', 'log', 'XTick', tauN, 'XTickLabel', {'25', '50', '100', '200'}); xlim(ax, [20 230]); ylim(ax, [1.4 3.5]);
xlabel(ax, '\tau_N (ms)', 'FontSize', S.fs); ylabel(ax, 'loss of quiet state, g_E^{tot}', 'FontSize', S.fs);
sa(ax, 'c  edge vs NMDA kinetics');
legend(ax, hc, {'f_N = 0.25', '0.5', '0.75'}, 'Location', 'east', 'Box', 'off');
% d
ax = subplot(2, 3, 4); hold(ax, 'on'); Mg = [0.5 1 2];
for s = 1:3
    e = arrayfun(@(m) edge(ld(s, 0.5, m), 3), Mg);
    plot(ax, Mg, e, ['-' S.mk{s}], 'Color', S.col(s,:), 'MarkerFaceColor', S.col(s,:), 'MarkerSize', S.ms, 'LineWidth', 1.6);
end
set(ax, 'XScale', 'log', 'XTick', Mg, 'XTickLabel', {'0.5', '1', '2'}); xlim(ax, [0.4 2.5]); ylim(ax, [1 3.5]);
xlabel(ax, '[Mg^{2+}] (mM)', 'FontSize', S.fs); ylabel(ax, 'loss of quiet state, g_E^{tot}', 'FontSize', S.fs);
sa(ax, 'd  edge vs Mg^{2+} (f_N = 0.5)');
legend(ax, {'connectivity 1', '2', '3'}, 'Location', 'northwest', 'Box', 'off');
% e
ax = subplot(2, 3, 5); hold(ax, 'on'); E2 = load_scan('exp2_kinetics');
for k = 1:4
    sel = E2.pt(:,2) == tauN(k) & E2.has_fp == 1 & E2.rateE < 4;
    plot(ax, E2.pt(sel,1) + 0.015*(k-2.5), E2.lam_aut(sel), S.mk{k}, 'Color', S.col(k,:), ...
         'MarkerFaceColor', S.col(k,:), 'MarkerSize', S.ms);
    line(ax, [1.1 1.6], -1000/tauN(k)*[1 1], 'Color', S.col(k,:), 'LineStyle', '--', 'LineWidth', 1.4);
end
xlim(ax, [1.1 1.6]); ylim(ax, [-45 0]);
xlabel(ax, 'g_E^{tot} (\times g_L)', 'FontSize', S.fs); ylabel(ax, '\lambda_{aut} (s^{-1})', 'FontSize', S.fs);
sa(ax, 'e  kinetics set the relaxation rate');
% f
ax = subplot(2, 3, 6); hold(ax, 'on');
E3 = load_scan('exp3_magnesium'); Mg4 = [0.25 0.5 1 2];
for k = 1:4
    sel = abs(E3.pt(:,2) - Mg4(k)) < 1e-9 & E3.has_fp == 1;
    hd(k) = plot(ax, [NaN; E3.rhoK_noNSC(sel)], [NaN; E3.rhoK(sel)], S.mk{k}, 'Color', S.col(k,:), ...
                 'MarkerFaceColor', S.col(k,:), 'MarkerSize', S.ms); %#ok<AGROW>
end
line(ax, [0.2 1.1], [0.2 1.1], 'Color', [0.3 0.3 0.3], 'LineStyle', '--', 'LineWidth', 1.1);
xlim(ax, [0.2 1.05]); ylim(ax, [0.2 1.05]);
xlabel(ax, '\rho_K without NSC', 'FontSize', S.fs); ylabel(ax, '\rho_K', 'FontSize', S.fs);
sa(ax, 'f  NMDA negative slope conductance');
legend(ax, hd, {'0.25 mM', '0.5', '1', '2'}, 'Location', 'southeast', 'Box', 'off');
save_fig(fh, 'Fig4_theory', 15.5, 10.5);
close(fh);
end

function e = pick_edge(L, t)
% first loss of stability of the quiet branch: Hopf if before the fold
if isfinite(L.gE_unstable(t)), e = L.gE_unstable(t);
elseif isfield(L, 'gE_end') && isfinite(L.gE_end), e = L.gE_end;
else, e = L.gEs(end); end   % branch not ended within the scanned range (lower bound)
end
