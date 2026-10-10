function figS3_restoration()
%FIGS3_RESTORATION  Restoration from actual states of the homeostatic cycle
%   (Supplementary Fig. S3; Experiment 11). Full input is restored at
%   epochs 41, 43, ..., 59 of the cycle, from the stored network state and
%   with the scaling factor theta in force at that moment.
%   Top: mean E rate after restoration (dashed: 4 Hz); bottom: memory
%   capacity for the restored input; columns: connectivities 1-3.
%   Conditions: AMPA only, NMDA (fN = 0.5), NMDA with 2 mM Mg2+ at
%   restoration, NMDA with 50% NMDA conductance block at restoration.
S = fig_style(); fh = figure('Visible', 'off', 'Color', 'w');
rd = fullfile(fileparts(mfilename('fullpath')), '..', 'results', 'exp11_restore');
yl = {'E rate after restoration (Hz)', 'MC after restoration'}; fl = {'rateE', 'MC'};
lab = {'AMPA only', 'NMDA', 'NMDA + 2 mM Mg^{2+}', 'NMDA + 50% block'};
for s = 1:3
    for p = 1:2
        ax = subplot(2, 3, (p - 1) * 3 + s); hold(ax, 'on'); hh = zeros(1, 4);
        for c = 1:2
            L = load(fullfile(rd, sprintf('seed%d_cond%d.mat', s, c))); nm = fieldnames(L.R);
            for q = 1:numel(nm)
                ci = c + q - 1; th = [L.R.(nm{q}).theta]; y = [L.R.(nm{q}).(fl{p})];
                dx = 0.03 * (ci - 2.5);
                hh(ci) = plot(ax, th + dx, y, S.mk{ci}, 'Color', S.col(ci,:), 'MarkerFaceColor', S.col(ci,:), 'MarkerSize', S.ms);
            end
        end
        if p == 1, refline_h(ax, 4); ylim(ax, [0 20]); else, ylim(ax, [0 5]); end
        xlim(ax, [0 4.2]); xlabel(ax, '\theta at restoration'); ylabel(ax, yl{p});
        style_axes(ax, sprintf('%s  connectivity %d', char('a' + (p - 1) * 3 + s - 1), s));
        if p == 1 && s == 1, legend(ax, hh, lab, 'Location', 'northwest', 'Box', 'off'); end
    end
end
save_fig(fh, 'FigS3_restoration', 14, 9.5);
close(fh);
end
