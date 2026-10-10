function exp2_kinetics(worker, nworkers)
%EXP2_KINETICS  NMDA decay time constant tau_N at fixed composition fN = 0.5.
%   Tests the prediction that tau_N does not move the static instability
%   (the static gain K has no time constant) but rescales the dynamics.
if nargin < 1, worker = 1; nworkers = 1; end
gE = 1.2:0.2:3.0; tauN = [25 50 100 200]; seeds = 1:2;
pts = [];
for s = seeds, for t = tauN, for g = gE, pts(end+1, :) = [g t s]; end, end, end %#ok<AGROW>
make_P = @(p) setfield(set_composition(base_params(), 0.5, p(1), -55), 'tau_N', p(2));
run_grid('exp2_kinetics', pts, make_P, scan_cfg(), worker, nworkers);
end
