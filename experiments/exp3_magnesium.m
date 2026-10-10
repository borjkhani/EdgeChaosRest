function exp3_magnesium(worker, nworkers)
%EXP3_MAGNESIUM  Extracellular Mg2+ at fixed NMDA and AMPA conductances.
%   Conductances are set by set_composition at the physiological [Mg2+] =
%   1 mM (fN = 0.5) and then [Mg2+] is changed, so that lowering Mg2+
%   unblocks NMDA receptors (raises B(V)) without changing gA or gN.
if nargin < 1, worker = 1; nworkers = 1; end
gE = 1.2:0.2:3.0; Mg = [0.25 0.5 1 2]; seeds = 1:2;
pts = [];
for s = seeds, for m = Mg, for g = gE, pts(end+1, :) = [g m s]; end, end, end %#ok<AGROW>
make_P = @(p) setfield(set_composition(base_params(), 0.5, p(1), -55), 'Mg', p(2));
run_grid('exp3_magnesium', pts, make_P, scan_cfg(), worker, nworkers);
end
