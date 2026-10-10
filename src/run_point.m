function R = run_point(P, seed, cfg)
%RUN_POINT  All measures for one parameter set and one connectivity seed.
%   Builds the network (connectivity seed = seed), evaluates dynamics and
%   memory (evaluate_point) and, if the noise-free autonomous network
%   settles on a fixed point, the static-gain criterion (gain_criterion).
if nargin < 3, cfg = struct(); end
C = build_connectivity(P, seed);
R = evaluate_point(P, C, seed, cfg);
[x, conv] = relax_to_fp(P, C, [], 250);
R.fp_conv = conv;
R.has_fp  = conv < 1e-4;
if R.has_fp
    S = gain_criterion(P, C, x);
    f = fieldnames(S);
    for q = 1:numel(f), R.(f{q}) = S.(f{q}); end
else
    R.rhoK = NaN; R.rhoK_noNSC = NaN; R.rhoK_noMg = NaN; R.maxReJ = NaN;
    R.leadComplex = false; R.G = NaN; R.NSC = NaN; R.V = NaN; R.dnu = NaN;
end
R.seed = seed;
end
