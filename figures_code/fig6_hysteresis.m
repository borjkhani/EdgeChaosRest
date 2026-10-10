function fig6_hysteresis()
%FIG6_HYSTERESIS  Steady-state rate vs homeostatic scaling theta (exp. 5).
%   Rows: intact feedforward drive (f_in = 1) and degenerated retina
%   (f_in = 0.25). Columns: connectivity seeds 1 and 2. Filled markers and
%   solid lines: start from a low state; open markers, dashed: start from a
%   high state (differences = bistability). Dashed horizontal line: mean E
%   rate of the healthy network (homeostatic set point).
S = fig_style(); fh = figure('Visible', 'off', 'Color', 'w');
here = fileparts(mfilename('fullpath'));
names = {'AMPA only', 'NMDA, 1 mM Mg^{2+}', 'NMDA, 2 mM Mg^{2+}', 'NMDA, 50% block'};
rowlab = {'intact input', 'degenerated input (f_{in} = 0.25)'};
p = 0;
for kf = 1:2
    for seed = 1:2
        p = p + 1; ax = subplot(2, 2, p); hold(ax, 'on');
        for c = 1:4
            L = load(fullfile(here, '..', 'results', 'exp5_hysteresis', sprintf('cond%d_seed%d.mat', c, seed)));
            dx = 0.012 * (c - 2.5);
            h(c) = plot(ax, L.theta + dx, L.rate(:, kf, 1), ['-' S.mk{c}], 'Color', S.col(c,:), ...
                'MarkerFaceColor', S.col(c,:), 'MarkerSize', S.ms, 'LineWidth', S.lw); %#ok<AGROW>
            plot(ax, L.theta + dx, L.rate(:, kf, 2), ['--' S.mk{c}], 'Color', S.col(c,:), ...
                'MarkerSize', S.ms, 'LineWidth', 0.8);
        end
        xlim(ax, [0.3 3.1]); ylim(ax, [0 20]); refline_h(ax, 2.13);
        xlabel(ax, 'scaling factor \theta'); ylabel(ax, 'mean E rate (Hz)');
        style_axes(ax, sprintf('%s  %s, seed %d', char('a' + p - 1), rowlab{kf}, seed));
        if p == 1, legend(ax, h, names, 'Location', 'northwest', 'Box', 'off'); end
    end
end
save_fig(fh, 'Fig6_hysteresis', 17, 12);
close(fh);
end
