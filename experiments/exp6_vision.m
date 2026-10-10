function exp6_vision(seed, cond_idx)
%EXP6_VISION  Retinal deafferentation -> homeostatic scaling -> restoration.
%
%   V1 receives feedforward (retino-geniculate) drive
%   f_in * (I_ff + sig_u * w * u(t)).
%     1. Healthy      : f_in = 1. Mean E rate = homeostatic set point.
%     2. Degeneration : f_in = 0.25. Population-level multiplicative scaling
%                       theta of all excitatory (AMPA+NMDA) synapses tries to
%                       restore the set point (homeostatic_scaling, 60 epochs).
%                       theta_dep = mean theta over the last 20 epochs.
%     3. Restoration  : f_in = 0.5 (partial) or 1 (full), theta frozen at
%                       theta_dep (restoration is fast relative to
%                       homeostasis). Measured: rate, memory capacity for the
%                       restored input, conditional Lyapunov exponent.
%     4. Interventions at restoration (theta frozen): [Mg2+] = 2 mM, or a
%                       50% block of NMDA conductance (memantine-like).
%     5. Re-scaling   : homeostasis continues under restored input (40 epochs).
%
%   cond_idx: 1 AMPA only, 2 NMDA (fN = 0.5) at 1 mM Mg2+ (condition_params).
%   Saves results/exp6_vision/seed<seed>_cond<cond_idx>.mat
if nargin < 1, seed = 1; end
if nargin < 2, cond_idx = 2; end
outdir = fullfile('results', 'exp6_vision');
if ~exist(outdir, 'dir'), mkdir(outdir); end
fn = fullfile(outdir, sprintf('seed%d_cond%d.mat', seed, cond_idx));

f_deg = 0.25;  f_rest = [0.5 1.0];
n_dep = 60; n_avg = 20; n_rec = 40;
hcfg = struct('mode', 'global', 'eta', 0.3, 'T_epoch', 200, 'seed', seed);
cfg = scan_cfg(); cfg.do_aut = false;

P = condition_params(cond_idx);
C = build_connectivity(P, seed);
one = ones(C.N, 1);

% ---- 1. healthy
randn('state', 500 + seed);
o = simulate_net(P, C, struct('u', randn(500, 1), 'noise', true, 'seed', 600 + seed));
r_target = P.r_max * mean(o.R(:, 101:end), 2);
Out.healthy = measure(P, C, one, seed, cfg);
Out.r_target = r_target;

% ---- 2. degeneration + homeostatic scaling
Pd = P; Pd.f_in = f_deg;
Out.deg_before = measure(Pd, C, one, seed, cfg);
[~, Hd] = homeostatic_scaling(Pd, C, one, r_target, n_dep, hcfg);
Out.hist_dep = Hd;
theta_dep = mean(Hd.theta(end-n_avg+1:end));
Out.theta_dep = theta_dep;
Out.deg_after = measure(Pd, C, theta_dep * one, seed, cfg);

% ---- 3-5. restoration
for k = 1:numel(f_rest)
    Pr = P; Pr.f_in = f_rest(k);
    Rk = struct();
    Rk.f_rest   = f_rest(k);
    Rk.frozen   = measure(Pr, C, theta_dep * one, seed, cfg);
    Pm = Pr; Pm.Mg = 2.0;
    Rk.high_mg  = measure(Pm, C, theta_dep * one, seed, cfg);
    Pb = Pr; Pb.gN = 0.5 * Pr.gN;
    Rk.nmda_block = measure(Pb, C, theta_dep * one, seed, cfg);
    [sc_rec, Hr] = homeostatic_scaling(Pr, C, theta_dep * one, r_target, n_rec, ...
                                       setfield(hcfg, 'seed', seed + 50));
    Rk.hist_rec = Hr;
    Rk.theta_rec = mean(Hr.theta(end-9:end));
    Rk.rescaled = measure(Pr, C, Rk.theta_rec * one, seed, cfg);
    Out.rest(k) = Rk;
end
Out.seed = seed; Out.cond = cond_idx; Out.f_deg = f_deg; Out.P = P;
save(fn, 'Out');
fprintf('saved %s (theta_dep = %.2f)\n', fn, theta_dep);
end

function M = measure(P, C, scaleE, seed, cfg)
%MEASURE  Rate, memory capacity and conditional Lyapunov exponent; the
%   network starts from a low-activity state.
cfg.scaleE = scaleE;
M = evaluate_point(P, C, seed, cfg);
end
