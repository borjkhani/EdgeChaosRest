function exp11_restore_states(seed, cond_idx)
%EXP11_RESTORE_STATES  Restoration from actual states of the homeostatic
%   cycle. The degeneration and scaling of Experiment 6 is repeated (same
%   random sequences, so the same trajectory) while the network state is
%   stored after every epoch. Full input (f_in = 1) is then restored at
%   epochs 41, 43, ..., 59, starting from the stored network state and with
%   the scaling factor in force at that moment (frozen during the test).
%   For the NMDA condition, 2 mM Mg2+ and 50% NMDA conductance block are
%   also applied at restoration.
%   Saves results/exp11_restore/seed<s>_cond<c>.mat
outdir = fullfile('results', 'exp11_restore');
if ~exist(outdir, 'dir'), mkdir(outdir); end
fn = fullfile(outdir, sprintf('seed%d_cond%d.mat', seed, cond_idx));
if exist(fn, 'file'), return; end
P = condition_params(cond_idx); C = build_connectivity(P, seed); one = ones(C.N, 1);
randn('state', 500 + seed);
o = simulate_net(P, C, struct('u', randn(500, 1), 'noise', true, 'seed', 600 + seed));
r_target = P.r_max * mean(o.R(:, 101:end), 2);
hcfg = struct('mode', 'global', 'eta', 0.3, 'T_epoch', 200, 'seed', seed, 'save_x', true);
Pd = P; Pd.f_in = 0.25;
[~, H] = homeostatic_scaling(Pd, C, one, r_target, 60, hcfg);
cfg = scan_cfg(); cfg.do_aut = false;
ep = 41:2:59; nm = {'full', 'high_mg', 'nmda_block'};
if cond_idx == 1, nm = nm(1); end
R = struct();
for q = 1:numel(nm)
    for k = 1:numel(ep)
        e = ep(k); Pr = P; Pr.f_in = 1;
        if q == 2, Pr.Mg = 2; end
        if q == 3, Pr.gN = 0.5 * Pr.gN; end
        cfg.scaleE = H.theta(e) * one; cfg.x0 = H.x{e};
        M = evaluate_point(Pr, C, seed, cfg);
        R.(nm{q})(k) = struct('epoch', e, 'theta', H.theta(e), 'rate_before', H.rateE(e), ...
            'rateE', M.rateE, 'MC', M.MC, 'lam_cond', M.lam_cond);
    end
end
H = rmfield(H, 'x');
save(fn, 'R', 'H', 'ep', 'seed', 'cond_idx');
fprintf('saved %s\n', fn);
end
