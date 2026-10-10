function exp5_hysteresis(cond_idx, seed)
%EXP5_HYSTERESIS  Bistability of the scaled network under intact and lost input.
%   For a grid of homeostatic scaling factors theta (all excitatory
%   synapses multiplied by theta), the network is started from a low
%   (V_init = -65 mV) and from a high (V_init = -35 mV) state, with
%   feedforward drive f_in = 1 (intact) or f_in = 0.25 (degenerated
%   retina). The steady-state mean E rate, its temporal fluctuation and the
%   conditional Lyapunov exponent are recorded. Differences between the two
%   branches reveal hysteresis (bistability) and the point where the low,
%   input-sensitive branch terminates.
%   Conditions: see condition_params.m.
if nargin < 2, seed = 1; end
outdir = fullfile('results', 'exp5_hysteresis');
if ~exist(outdir, 'dir'), mkdir(outdir); end
fn = fullfile(outdir, sprintf('cond%d_seed%d.mat', cond_idx, seed));
P = condition_params(cond_idx);
C = build_connectivity(P, seed);
theta = 0.4:0.2:3.0; fin = [1 0.25]; Vinit = [-65 -35];
rate = nan(numel(theta), 2, 2); sdr = rate; lam = rate;
for it = 1:numel(theta)
    for kf = 1:2
        for kv = 1:2
            Pk = P; Pk.f_in = fin(kf);
            randn('state', 300 + it);
            o = simulate_net(Pk, C, struct('u', randn(200, 1), 'noise', true, 'lyap', true, ...
                'lyap_skip', 100, 'seed', seed * 31 + it, 'scaleE', theta(it) * ones(C.N, 1), ...
                'V_init', Vinit(kv)));
            rE = Pk.r_max * mean(o.R(C.isE, 101:end), 1);
            rate(it, kf, kv) = mean(rE); sdr(it, kf, kv) = std(rE); lam(it, kf, kv) = o.lam;
        end
    end
    fprintf('cond%d s%d theta=%.1f  intact: %5.2f/%5.2f   deg: %5.2f/%5.2f Hz\n', cond_idx, seed, theta(it), ...
        rate(it,1,1), rate(it,1,2), rate(it,2,1), rate(it,2,2));
end
save(fn, 'theta', 'fin', 'Vinit', 'rate', 'sdr', 'lam', 'cond_idx', 'seed', 'P');
end
