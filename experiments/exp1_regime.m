function exp1_regime(worker, nworkers)
%EXP1_REGIME  Excitatory gain x AMPA/NMDA composition (Figures 2-4).
%   gE_tot : total excitatory conductance at V_ref (x g_L)
%   fN     : NMDA fraction of that conductance (set_composition)
%   10 connectivity seeds per point (seeds 1-3 are also used in Exp. 2-3b).
if nargin < 1, worker = 1; nworkers = 1; end
gE = 1.2:0.2:3.0; fN = [0 0.25 0.5 0.75]; seeds = 1:10;
pts = [];
for s = seeds, for f = fN, for g = gE, pts(end+1, :) = [g f s]; end, end, end %#ok<AGROW>
make_P = @(p) set_composition(base_params(), p(2), p(1), -55);
run_grid('exp1_regime', pts, make_P, scan_cfg(), worker, nworkers);
end
