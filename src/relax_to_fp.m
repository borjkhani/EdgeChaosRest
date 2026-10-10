function [x, conv, R] = relax_to_fp(P, C, scaleE, T_bins, x0)
%RELAX_TO_FP  Run the noise-free network with constant input (u = 0) and
%   return the final state. conv is the max |dV| over the last bin (mV),
%   which is small when the trajectory has settled on a fixed point.
if nargin < 3, scaleE = []; end
if nargin < 4 || isempty(T_bins), T_bins = 150; end
opt = struct('T_bins', T_bins, 'noise', false, 'scaleE', scaleE, 'seed', 7);
if nargin > 4 && ~isempty(x0), opt.x0 = x0; end
o = simulate_net(P, C, opt);
x = o.x; R = o.R;
conv = max(abs(P.theta + P.kappa * log(R(:, end) ./ (1 - R(:, end))) ...
             - P.theta - P.kappa * log(R(:, end-1) ./ (1 - R(:, end-1)))));
end
