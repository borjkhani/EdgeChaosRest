function exp8_homeostasis_robustness(seed, cond_idx, variant)
%EXP8_HOMEOSTASIS_ROBUSTNESS  Does homeostatic cycling after degeneration
%   depend on the update gain, the step size or the bounds of the scaling
%   rule? Degeneration (f_in = 0.25) followed by population-level scaling,
%   as in Experiment 6, with the rule changed:
%     'eta0.3'  : reference (eta = 0.3, 4-s epochs, 60 epochs), repeated
%     'eta0.1'  : eta = 0.1, 120 epochs
%     'eta0.03' : eta = 0.03, 200 epochs (250 for connectivity 1)
%     'tanh'    : bounded step, theta changes by at most exp(+-0.1) per epoch,
%                 120 epochs
%     'cont'    : quasi-continuous adaptation: 1-s epochs, eta = 0.0075
%                 (the same adaptation rate per second as eta = 0.03),
%                 700 epochs
%     'floor'   : eta = 0.3 with theta_min = 0.5 instead of 0.25, 60 epochs
%   The healthy set point is measured exactly as in Experiment 6.
%   Saves results/exp8_homeostasis/seed<s>_cond<c>_<variant>.mat
outdir = fullfile('results', 'exp8_homeostasis');
if ~exist(outdir, 'dir'), mkdir(outdir); end
fn = fullfile(outdir, sprintf('seed%d_cond%d_%s.mat', seed, cond_idx, variant));
if exist(fn, 'file'), return; end
P = condition_params(cond_idx); C = build_connectivity(P, seed); one = ones(C.N, 1);
randn('state', 500 + seed);
o = simulate_net(P, C, struct('u', randn(500, 1), 'noise', true, 'seed', 600 + seed));
r_target = P.r_max * mean(o.R(:, 101:end), 2);
h = struct('mode', 'global', 'eta', 0.3, 'T_epoch', 200, 'seed', seed);
switch variant
    case 'eta0.3',  n = 60;
    case 'eta0.1',  h.eta = 0.1;  n = 120;
    case 'eta0.03', h.eta = 0.03; n = 200;
    case 'tanh',    h.eta = 0.1;  h.rule = 'tanh'; n = 120;
    case 'cont',    h.eta = 0.0075; h.T_epoch = 50; n = 700;
    case 'floor',   h.theta_min = 0.5; n = 60;
end
Pd = P; Pd.f_in = 0.25;
t0 = tic;
[~, H] = homeostatic_scaling(Pd, C, one, r_target, n, h);
target = mean(r_target(C.isE));
save(fn, 'H', 'target', 'variant', 'h', 'seed', 'cond_idx');
fprintf('saved %s (%.0f s)\n', fn, toc(t0));
end
