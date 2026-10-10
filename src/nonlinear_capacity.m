function [NLC, parts] = nonlinear_capacity(R, u, n_wash, n_train, K)
%NONLINEAR_CAPACITY  Quadratic information-processing capacity
%   (Dambre et al. 2012) from sampled states, with the same ridge readout,
%   training and test split as memory_capacity.m.
%   Targets (all degree-2 functions of the past input up to lag K):
%     z_k   = (u(t-k)^2 - 1)/sqrt(2)          k = 1..K        (Legendre P2)
%     z_kl  = u(t-k) u(t-l)                   1 <= k < l <= K
%   NLC = sum over targets of the squared test correlation.
%   parts.square = capacity of the z_k, parts.cross = capacity of the z_kl.
u = u(:); T = size(R, 2);
X = [R; ones(1, T)]';
tr = (n_wash + K + 1):(n_wash + n_train);
te = (n_wash + n_train + 1):T;
lam = 1e-4 * mean(var(R(:, tr), 0, 2)) * numel(tr);
Xtr = X(tr, :); Mi = (Xtr' * Xtr + lam * eye(size(X, 2))) \ Xtr';   % reused for every target
Xte = X(te, :);
cap = @(z) corr2(Xte * (Mi * z(tr)), z(te));
sq = 0; cr = 0;
for k = 1:K
    z = zeros(T, 1); z(k+1:end) = (u(1:end-k).^2 - 1) / sqrt(2);
    sq = sq + cap(z);
    for l = k+1:K
        z = zeros(T, 1); z(l+1:end) = u(l-k+1:end-k) .* u(1:end-l);
        cr = cr + cap(z);
    end
end
NLC = sq + cr; parts = struct('square', sq, 'cross', cr);
end

function c = corr2(a, b)
c = corrcoef(a, b);
if all(isfinite(c(:))), c = max(0, c(1, 2))^2; else, c = 0; end
end
