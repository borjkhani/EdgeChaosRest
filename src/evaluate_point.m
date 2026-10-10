function M = evaluate_point(P, C, seed, cfg)
%EVALUATE_POINT  Dynamical and computational measures for one network.
%
%   M = evaluate_point(P, C, seed, cfg)
%
%   Runs (1) an autonomous, noise-free simulation (u = 0) to obtain the
%   autonomous maximal Lyapunov exponent lam_aut, and (2) an input-driven,
%   noisy simulation with i.i.d. Gaussian input u(t) ~ N(0,1) held for
%   P.bin_ms, from which the conditional Lyapunov exponent lam_cond
%   (exponent of the driven system, same input and noise realisation) and
%   the linear memory capacity MC are computed.
%
%   cfg fields (defaults in brackets): T_aut [400], skip_aut [150],
%   n_wash [200], n_train [2000], n_test [1000], K [60], scaleE [],
%   do_aut [true], x0 [] (initial state of the driven run), extra [false].
%   With extra = true the driven run also yields the quadratic (nonlinear)
%   memory capacity NLC (nonlinear_capacity.m) and, as chance baselines, the
%   linear and quadratic capacities for an independent input sequence that
%   was never presented (MC_shuf, NLC_shuf).
if nargin < 4, cfg = struct(); end
d = struct('T_aut', 400, 'skip_aut', 150, 'n_wash', 200, 'n_train', 2000, ...
           'n_test', 1000, 'K', 60, 'scaleE', [], 'do_aut', true, 'x0', [], 'extra', false);
f = fieldnames(d);
for q = 1:numel(f), if ~isfield(cfg, f{q}), cfg.(f{q}) = d.(f{q}); end, end

M = struct();
if cfg.do_aut
    oa = simulate_net(P, C, struct('T_bins', cfg.T_aut, 'lyap', true, ...
        'lyap_skip', cfg.skip_aut, 'seed', seed, 'scaleE', cfg.scaleE));
    M.lam_aut  = oa.lam;
    Ra = oa.R(:, cfg.skip_aut+1:end);
    M.sd_aut   = mean(std(Ra, 0, 2));             % temporal fluctuation per unit
    M.rateE_aut = oa.rateE;
else
    M.lam_aut = NaN; M.sd_aut = NaN; M.rateE_aut = NaN;
end

T = cfg.n_wash + cfg.n_train + cfg.n_test;
randn('state', 1000 + seed);
u = randn(T, 1);
od = simulate_net(P, C, struct('u', u, 'noise', true, 'lyap', true, ...
    'lyap_skip', cfg.n_wash, 'seed', 2000 + seed, 'scaleE', cfg.scaleE, 'x0', cfg.x0));
[M.MC, M.MCk] = memory_capacity(od.R, u, cfg.n_wash, cfg.n_train, cfg.K);
if cfg.extra
    randn('state', 3000 + seed); us = randn(T, 1);    % never presented
    M.MC_shuf = memory_capacity(od.R, us, cfg.n_wash, cfg.n_train, cfg.K);
    M.NLC = nonlinear_capacity(od.R, u, cfg.n_wash, cfg.n_train, 10);
    M.NLC_shuf = nonlinear_capacity(od.R, us, cfg.n_wash, cfg.n_train, 10);
end
M.lam_cond = od.lam;
M.rateE = od.rateE; M.rateI = od.rateI;
Rd = od.R(:, cfg.n_wash+1:end);
pm = mean(Rd, 1);
M.sync = var(pm) / mean(var(Rd, 0, 2));            % population synchrony
M.sd_drv = mean(std(Rd, 0, 2));
% memory time: lag at which MC_k first drops below 0.5 (in ms)
k50 = find(M.MCk < 0.5, 1);
if isempty(k50), k50 = cfg.K + 1; end
M.t_mem50 = (k50 - 1) * P.bin_ms;
end
