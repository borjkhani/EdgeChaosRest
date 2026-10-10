function exp10_nonlinear(worker, nworkers)
%EXP10_NONLINEAR  Memory controls for every run of Experiment 1: quadratic
%   (nonlinear) memory capacity and chance baselines for an input sequence
%   that was never presented (evaluate_point with extra = true). The driven
%   simulation is identical to Experiment 1 (same seeds), so MC is the same.
if nargin < 1, worker = 1; nworkers = 1; end
gE = 1.2:0.2:3.0; fN = [0 0.25 0.5 0.75]; seeds = 1:10;
pts = [];
for s = seeds, for f = fN, for g = gE, pts(end+1, :) = [g f s]; end, end, end %#ok<AGROW>
make_P = @(p) set_composition(base_params(), p(2), p(1), -55);
cfg = scan_cfg(); cfg.do_aut = false; cfg.extra = true;
run_grid('exp10_nonlinear', pts, make_P, cfg, worker, nworkers);
end
