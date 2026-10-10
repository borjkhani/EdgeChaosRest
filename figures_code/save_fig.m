function save_fig(fh, name, w_cm, h_cm)
%SAVE_FIG  Export figure fh to ../figures/<name>.pdf and .png (300 dpi).
if nargin < 3, w_cm = 17; end
if nargin < 4, h_cm = 12; end
outdir = fullfile(fileparts(mfilename('fullpath')), '..', 'figures');
if ~exist(outdir, 'dir'), mkdir(outdir); end
k = 1.6;   % render on a larger page so that 9-pt text is proportionate after LaTeX scaling
set(fh, 'PaperUnits', 'centimeters', 'PaperSize', k * [w_cm h_cm], ...
        'PaperPosition', [0 0 k * w_cm k * h_cm]);
print(fh, '-dpdf', fullfile(outdir, [name '.pdf']));
print(fh, '-dpng', '-r300', fullfile(outdir, [name '.png']));
end
