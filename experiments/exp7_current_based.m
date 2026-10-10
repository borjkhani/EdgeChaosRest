function exp7_current_based(worker, nworkers)
%EXP7_CURRENT_BASED  Robustness: AMPA and GABA_A current-based (fixed driving
%   force at V_ref), NMDA with full Mg2+ voltage dependence. A global
%   coupling J scales excitation (gE_tot = J) and inhibition (gG = 2.2 J).
%   Denser, less heterogeneous wiring (KE = 40, KI = 10, CV = 1) and a
%   softer rate function (kappa = 4 mV) keep units in their graded range.
if nargin < 1, worker = 1; nworkers = 1; end
pts = [];
for s = 1:2, for f = [0 0.5], for J = 0.6:0.4:2.6, pts(end+1, :) = [J f s]; end, end, end %#ok<AGROW>
run_grid('exp7_current_based', pts, @make_P, scan_cfg(), worker, nworkers);
end

function P = make_P(p)
P = default_params(); P.syn_mode = 'current'; P.kappa = 4; P.KE = 40; P.KI = 10;
P.cv_w = 1; P.tau_G = 5; P.I_bg = 7; P.I_ff = 6; P.sig_u = 1.5; P.sig_eta = 0.3; P.V_ref = -55;
P = set_composition(P, p(2), p(1), -55);
P.gG = 2.2 * p(1);
end
