function exp9b_replicates()
%EXP9B_REPLICATES  Separates the effect of the integration step from that of
%   the input and noise realisation at the points of exp9_numerics: each
%   point is re-run with 4 further input/noise realisations at the standard
%   step (0.5 ms) and at half the step (0.25 ms). Recorded: mean E rate,
%   memory capacity and conditional exponent.
%   Saves results/exp9_numerics/replicates.mat
if exist(fullfile('results', 'exp9_numerics', 'replicates.mat'), 'file'), return; end
L = load(fullfile('results', 'exp9_numerics', 'numerics.mat'));
cfg = scan_cfg(); cfg.do_aut = false; rep = struct([]);
for r = 1:numel(L.res)
    p = L.res(r).pt; P = set_composition(base_params(), p(2), p(1), -55); C = build_connectivity(P, p(3));
    for h = 1:2
        Ph = P; if h == 2, Ph.dt = 0.25; end
        for k = 1:4
            M = evaluate_point(Ph, C, p(3) + 100 * k, cfg);
            rep(r).rate(h, k) = M.rateE; rep(r).MC(h, k) = M.MC; rep(r).lam_cond(h, k) = M.lam_cond;
        end
    end
    rep(r).pt = p;
    fprintf('%s: active fraction dt 0.5 %.2f, dt 0.25 %.2f | MC %.2f vs %.2f | lamC %.2f vs %.2f\n', mat2str(p, 3), ...
        mean(rep(r).rate(1,:) >= 4), mean(rep(r).rate(2,:) >= 4), mean(rep(r).MC(1,:)), mean(rep(r).MC(2,:)), mean(rep(r).lam_cond(1,:)), mean(rep(r).lam_cond(2,:)));
end
save(fullfile('results', 'exp9_numerics', 'replicates.mat'), 'rep');
end
