function exp4_tau_heterogeneity(worker, nworkers)
%EXP4_TAU_HETEROGENEITY  Homogeneous vs heterogeneous NMDA decay constants.
%   A GluN2A/GluN2B mixture is modelled by drawing each presynaptic unit's
%   tau_N log-uniformly in [30, 330] ms (mean ~125 ms) and compared with a
%   homogeneous tau_N = 100 ms, at fN = 0.5 and 0.75, two couplings
%   (gE_tot = 1.2, 1.4) and three connectivity seeds.
if nargin < 1, worker = 1; nworkers = 1; end
pts = [];
for s = 1:3, for f = [0.5 0.75], for g = [1.2 1.4], for het = [0 1]
    pts(end+1, :) = [g f het s]; %#ok<AGROW>
end, end, end, end
run_grid('exp4_tau_heterogeneity', pts, @make_P, scan_cfg(), worker, nworkers);
end

function P = make_P(p)
P = set_composition(base_params(), p(2), p(1), -55);
if p(3) == 1
    rand('state', 99 + p(4));
    P.tau_N = exp(log(30) + rand(P.NE, 1) * (log(330) - log(30)));
end
end
