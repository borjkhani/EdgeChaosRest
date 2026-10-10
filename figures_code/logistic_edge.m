function [g50, ci, b, pfit] = logistic_edge(x, y, seed, nboot, xgrid)
%LOGISTIC_EDGE  Coupling at which the active state becomes more likely than
%   the quiet one, from a logistic fit  logit P(active) = b(1) + b(2) x.
%   x     : coupling of each run
%   y     : 1 if the run is in the active state, 0 if quiet
%   seed  : connectivity of each run; the 95% confidence interval is
%           obtained by resampling connectivities (cluster bootstrap), because
%           runs of the same connectivity are not independent
%   nboot : bootstrap resamples (default 2000)
%   xgrid : optional grid on which the fitted probability is returned
%   g50   : -b(1)/b(2), the coupling with P(active) = 0.5
%   A weak ridge penalty on the slope keeps the fit finite when a resample
%   happens to separate quiet and active runs perfectly.
if nargin < 4 || isempty(nboot), nboot = 2000; end
if nargin < 5, xgrid = []; end
x = x(:); y = double(y(:)); seed = seed(:);
b = fit(x, y); g50 = -b(1) / b(2);
pfit = 1 ./ (1 + exp(-(b(1) + b(2) * xgrid(:))));
us = unique(seed); gb = nan(nboot, 1);
st = rand('state'); rand('state', 1);   % reproducible bootstrap
for r = 1:nboot
    pick = us(ceil(rand(numel(us), 1) * numel(us)));
    idx = cell2mat(arrayfun(@(s) find(seed == s), pick, 'UniformOutput', false));
    bb = fit(x(idx), y(idx)); gb(r) = -bb(1) / bb(2);
end
rand('state', st);
gb = sort(gb(isfinite(gb)));
ci = [gb(max(1, round(0.025 * numel(gb)))), gb(round(0.975 * numel(gb)))];
end

function b = fit(x, y)
% penalised IRLS (Newton) for logistic regression
X = [ones(size(x)) x]; b = [0; 0]; L = diag([0 1e-2]);
for it = 1:100
    p = 1 ./ (1 + exp(-X * b)); W = p .* (1 - p) + 1e-9;
    g = X' * (y - p) - L * b; H = X' * (X .* W) + L;
    db = H \ g; b = b + db;
    if max(abs(db)) < 1e-10, break; end
end
end
