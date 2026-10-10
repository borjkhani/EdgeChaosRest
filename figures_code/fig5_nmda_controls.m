function fig5_nmda_controls()
%FIG5_NMDA_CONTROLS  Mg2+ moves the edge of the quiet state; NMDA kinetics
%   change exponents but not memory.
%   a,b  mean E rate and memory capacity vs gE_tot for four [Mg2+]
%        (experiment 3; filled = quiet state, open = active state).
%   c    memory capacity vs gE_tot for four NMDA decay constants
%        (experiment 2), same encoding.
%   d    conditional Lyapunov exponent vs gE_tot for the same runs.
%   e    linearised quiet state (Experiment 13): relative power of the slowest mode's
%        component of the input impulse response of the rates (lags
%        20-1000 ms), relative to the total, versus tau_N, for the two
%        connectivities (markers) and two couplings (colours).
%   f    impulse-response power of the rates versus lag (connectivity 2,
%        gE_tot = 1.4) for tau_N = 25 and 200 ms; solid: full, dashed:
%        with the input to the slowest mode removed.
S = fig_style(); S.fs = 10; fh = figure('Visible', 'off', 'Color', 'w');
E3 = load_scan('exp3_magnesium'); E2 = load_scan('exp2_kinetics');
Mg = [0.25 0.5 1 2]; tauN = [25 50 100 200];
ax = subplot(2, 3, 1); hold(ax, 'on');
h = plot_states(ax, E3, 1, 2, Mg, 'rateE', S); refline_h(ax, 4);
style_axes(ax, 'a  [Mg^{2+}]: mean E rate (Hz)', S.fs);
ylim(ax, [0 44]); legend(ax, h, {'[Mg^{2+}] 0.25 mM', '0.5 mM', '1 mM', '2 mM'}, 'Location', 'northwest', 'Box', 'off');
ax2 = subplot(2, 3, 2); hold(ax2, 'on');
plot_states(ax2, E3, 1, 2, Mg, 'MC', S); style_axes(ax2, 'b  [Mg^{2+}]: memory capacity', S.fs);
ax3 = subplot(2, 3, 3); hold(ax3, 'on');
h3 = plot_states(ax3, E2, 1, 2, tauN, 'MC', S); style_axes(ax3, 'c  \tau_N: memory capacity', S.fs);
legend(ax3, h3, {'\tau_N = 25 ms', '50 ms', '100 ms', '200 ms'}, 'Location', 'northeast', 'Box', 'off');
ax4 = subplot(2, 3, 4); hold(ax4, 'on');
plot_states(ax4, E2, 1, 2, tauN, 'lam_cond', S); ylim(ax4, [-30 12]); refline_h(ax4, 0);
style_axes(ax4, 'd  \tau_N: \lambda_{cond} (s^{-1})', S.fs);
for a = [ax ax2 ax3 ax4], xlim(a, [1.1 3.1]); xlabel(a, 'g_E^{tot} (\times g_L)'); end
ylim(ax2, [0 5]); ylim(ax3, [0 7.5]);
% e, f: slowest mode and the input signal (Experiment 13)
L = load(fullfile(fileparts(mfilename('fullpath')), '..', 'results', 'exp13_slow_mode', 'slow_mode.mat'));
o = L.out; ax5 = subplot(2, 3, 5); hold(ax5, 'on'); he = [];
for sd = 1:2, for gi = 1:2
    g = [1.2 1.4]; k = o(:,1) == sd & abs(o(:,2) - g(gi)) < 1e-9;
    he(end+1) = plot(ax5, o(k,3), o(k,5), ['-' S.mk{sd}], 'Color', S.col(gi,:), 'MarkerFaceColor', S.col(gi,:), 'MarkerSize', S.ms, 'LineWidth', 1.2); %#ok<AGROW>
end, end
set(ax5, 'XScale', 'log', 'YScale', 'log', 'XTick', tauN, 'XTickLabel', {'25', '50', '100', '200'}); xlim(ax5, [20 230]); ylim(ax5, [5e-4 1.5]);
xlabel(ax5, '\tau_N (ms)'); ylabel(ax5, 'relative component power');
style_axes(ax5, 'e  slow-mode component power', S.fs);
legend(ax5, he, {'conn. 1, g_E^{tot} = 1.2', 'conn. 1, 1.4', 'conn. 2, 1.2', 'conn. 2, 1.4'}, 'Location', 'southwest', 'Box', 'off');
ax6 = subplot(2, 3, 6); hold(ax6, 'on'); hf = [];
for t = [25 200]
    k = find(o(:,1) == 2 & abs(o(:,2) - 1.4) < 1e-9 & o(:,3) == t); c = L.curves{k}; ci = 1 + (t == 200) * 3;
    hf(end+1) = plot(ax6, 1000 * L.lags, c(1,:), '-', 'Color', S.col(ci,:), 'LineWidth', S.lw); %#ok<AGROW>
    plot(ax6, 1000 * L.lags, c(2,:), '--', 'Color', S.col(ci,:), 'LineWidth', S.lw);
end
hf(end+1) = plot(ax6, NaN, NaN, '--', 'Color', S.ink2, 'LineWidth', S.lw);
set(ax6, 'YScale', 'log'); xlim(ax6, [0 400]); ylim(ax6, [1e-3 300]);
xlabel(ax6, 'lag (ms)'); ylabel(ax6, 'signal power ||h||^2');
style_axes(ax6, 'f  impulse response (conn. 2, 1.4)', S.fs);
legend(ax6, hf, {'\tau_N = 25 ms', '\tau_N = 200 ms', 'slow mode removed'}, 'Location', 'northeast', 'Box', 'off');
save_fig(fh, 'Fig5_nmda_controls', 15.5, 10.5);
close(fh);
end
