function [MC, MCk] = memory_capacity(R, u, n_wash, n_train, K, reg)
%MEMORY_CAPACITY  Linear memory capacity (Jaeger 2002) from sampled states.
%   R : N x T states (one column per input bin), u : T x 1 i.i.d. input.
%   Bins 1..n_wash are discarded; the next n_train bins train a ridge
%   readout for each lag k = 1..K, the remaining bins are used for testing.
%   MC_k = squared Pearson correlation between u(t-k) and its reconstruction
%   on the test set; MC = sum_k MC_k. A constant regressor is included.
%   reg is the ridge parameter relative to the mean state variance.
if nargin < 6, reg = 1e-4; end
u = u(:);
T = size(R, 2);
X = [R; ones(1, T)]';
MCk = zeros(K, 1);
tr  = (n_wash + 1):(n_wash + n_train);
te  = (n_wash + n_train + 1):T;
lam = reg * mean(var(R(:, tr), 0, 2)) * numel(tr);
for k = 1:K
    trk = tr(tr - k >= 1); tek = te;
    Xtr = X(trk, :);  ytr = u(trk - k);
    M   = Xtr' * Xtr + lam * eye(size(X, 2));
    wout = M \ (Xtr' * ytr);
    yhat = X(tek, :) * wout;  y = u(tek - k);
    c = corrcoef(yhat, y);
    if all(isfinite(c(:)))
        MCk(k) = max(0, c(1, 2))^2;
    end
end
MC = sum(MCk);
end
