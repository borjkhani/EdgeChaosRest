function [scaleE, hist] = homeostatic_scaling(P, C, scaleE, r_target, n_epochs, cfg)
%HOMEOSTATIC_SCALING  Cell-autonomous multiplicative synaptic scaling.
%
%   Excitatory (AMPA and NMDA) conductances onto unit i are multiplied by
%   theta_i = scaleE(i). After each epoch of driven activity the factor is
%   updated towards the unit's own target rate (Turrigiano-type scaling):
%
%       theta_i <- theta_i * exp( eta * (r*_i - r_i) / max(r*_i, r_floor) )
%
%   clipped to [theta_min, theta_max]. One epoch stands for the slow
%   (hours-days) homeostatic time scale; within an epoch theta is fixed.
%   Scaling acts on E and I units alike (both receive excitatory synapses
%   that scale), following the convention that all excitatory synapses
%   onto a cell scale together (AMPA and NMDA components alike).
%
%   With cfg.mode = 'global' a single factor theta (shared by all units)
%   is driven by the population-mean excitatory rate instead:
%       theta <- theta * exp( eta * (mean r*_E - mean r_E) / mean r*_E ).
%   This population-level version is used in the manuscript; the
%   cell-autonomous version is provided for comparison.
%
%   With cfg.rule = 'tanh' the relative error is passed through tanh before
%   the update, so that one epoch can change theta by at most a factor
%   exp(+-eta) (bounded-step variant used as a robustness check).
%
%   r_target : N x 1 target rates (Hz) measured in the healthy network
%   cfg      : mode ['global'], rule ['exp'], eta [0.3], r_floor [1 Hz],
%              theta_min [0.25], theta_max [6], T_epoch [200 bins], seed [1],
%              save_x [false] (store the network state after every epoch),
%              x0 [] (initial network state)
%   hist     : per-epoch mean theta, mean E rate, conditional Lyapunov
%              exponent and fraction of units within 20% of target.
if nargin < 6, cfg = struct(); end
d = struct('mode', 'global', 'rule', 'exp', 'eta', 0.3, 'r_floor', 1, 'theta_min', 0.25, ...
           'theta_max', 6, 'T_epoch', 200, 'seed', 1, 'save_x', false, 'x0', []);
f = fieldnames(d);
for q = 1:numel(f), if ~isfield(cfg, f{q}), cfg.(f{q}) = d.(f{q}); end, end
if numel(scaleE) == 1, scaleE = scaleE * ones(C.N, 1); end
hist = struct('theta', zeros(n_epochs, 1), 'rateE', zeros(n_epochs, 1), ...
              'lam', zeros(n_epochs, 1), 'frac_ok', zeros(n_epochs, 1));
x0 = cfg.x0;
if cfg.save_x, hist.x = cell(n_epochs, 1); end
for e = 1:n_epochs
    randn('state', cfg.seed * 7919 + e);
    u = randn(cfg.T_epoch, 1);
    sk = min(50, floor(cfg.T_epoch / 4));      % transient discarded within an epoch
    opt = struct('u', u, 'noise', true, 'lyap', true, 'lyap_skip', sk, ...
                 'seed', cfg.seed * 104729 + e, 'scaleE', scaleE, 'x0', x0);
    o = simulate_net(P, C, opt);
    x0 = o.x;
    r = P.r_max * mean(o.R(:, sk+1:end), 2);
    if strcmp(cfg.mode, 'global')
        rt = mean(r_target(C.isE));
        err = (rt - mean(r(C.isE))) / max(rt, cfg.r_floor);
    else
        err = (r_target - r) ./ max(r_target, cfg.r_floor);
    end
    if strcmp(cfg.rule, 'tanh'), err = tanh(err); end
    scaleE = min(max(scaleE .* exp(cfg.eta * err), cfg.theta_min), cfg.theta_max);
    if cfg.save_x, hist.x{e} = x0; end
    hist.theta(e)   = mean(scaleE);
    hist.rateE(e)   = mean(r(C.isE));
    hist.lam(e)     = o.lam;
    hist.frac_ok(e) = mean(abs(r - r_target) <= 0.2 * max(r_target, cfg.r_floor));
end
end
