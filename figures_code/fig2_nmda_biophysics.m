function fig2_nmda_biophysics()
%FIG2_NMDA_BIOPHYSICS  Mg2+ block, NMDA I-V relation and NMDA slope factor.
%   (a) B(V) for [Mg2+] = 0.5, 1, 2 mM (Jahr & Stevens 1990)
%   (b) normalised NMDA current B(V)(E_E - V)
%   (c) slope factor d/dV[B(V)(V - E_E)] = B + B'(V - E_E); negative values
%       are the negative slope conductance that lowers the local leak G_i.
addpath(fullfile(fileparts(mfilename('fullpath')), '..', 'src'));
S = fig_style(); P = default_params();
V = linspace(-90, 20, 400); Mg = [0.5 1 2];
fh = figure('Visible', 'off', 'Color', 'w');
for k = 1:3
    ax(k) = subplot(1, 3, k); hold(ax(k), 'on'); %#ok<AGROW>
end
for m = 1:numel(Mg)
    P.Mg = Mg(m); [B, dB] = mg_block(V, P);
    I = B .* (P.E_E - V);
    slope = B + dB .* (V - P.E_E);
    sel = 1:40:numel(V);
    plot(ax(1), V, B, '-', 'Color', S.col(m,:), 'LineWidth', S.lw);
    plot(ax(1), V(sel), B(sel), S.mk{m}, 'Color', S.col(m,:), 'MarkerFaceColor', S.col(m,:), 'MarkerSize', S.ms);
    plot(ax(2), V, I / max(I), '-', 'Color', S.col(m,:), 'LineWidth', S.lw);
    plot(ax(2), V(sel), I(sel) / max(I), S.mk{m}, 'Color', S.col(m,:), 'MarkerFaceColor', S.col(m,:), 'MarkerSize', S.ms);
    hs(m) = plot(ax(3), V, slope, '-', 'Color', S.col(m,:), 'LineWidth', S.lw); %#ok<AGROW>
    plot(ax(3), V(sel), slope(sel), S.mk{m}, 'Color', S.col(m,:), 'MarkerFaceColor', S.col(m,:), 'MarkerSize', S.ms);
end
xlabel(ax(1), 'V (mV)'); ylabel(ax(1), 'unblocked fraction B(V)');
xlabel(ax(2), 'V (mV)'); ylabel(ax(2), 'NMDA current (normalised)');
xlabel(ax(3), 'V (mV)'); ylabel(ax(3), 'B + B''(V - E_E)');
style_axes(ax(1), 'a  Mg^{2+} unblock'); style_axes(ax(2), 'b  NMDA I-V');
style_axes(ax(3), 'c  NMDA slope conductance');
for k = 1:3, xlim(ax(k), [-90 20]); end
ylim(ax(2), [-0.6 1.1]); ylim(ax(3), [-0.5 1.1]);
refline_h(ax(2), 0); refline_h(ax(3), 0);
legend(ax(1), hs, {'[Mg^{2+}] 0.5 mM', '1 mM', '2 mM'}, 'Location', 'northwest', 'Box', 'off');
save_fig(fh, 'Fig2_nmda_biophysics', 17, 6.2);
close(fh);
end
