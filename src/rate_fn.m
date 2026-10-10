function [nu, dnu] = rate_fn(V, P)
%RATE_FN  Normalised firing rate nu = phi(V) in (0,1) and its derivative.
nu  = 1 ./ (1 + exp(-(V - P.theta) ./ P.kappa));
dnu = nu .* (1 - nu) ./ P.kappa;
end
