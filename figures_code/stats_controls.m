function stats_controls()
%STATS_CONTROLS  Print every number quoted for the control analyses
%   (Section 4.6 and Supplementary Figs S1-S4): nonlinear memory and chance
%   baselines (Exp. 10), distance to the linear edge, numerical checks
%   (Exp. 9), arclength continuation (Exp. 3c), coexistence (Exp. 12), slow
%   mode (Exp. 13), homeostatic rules (Exp. 8) and restoration from actual
%   states (Exp. 11). Confidence intervals resample connectivities.
here = fileparts(mfilename('fullpath')); addpath(here); addpath(fullfile(here, '..', 'src'));
root = fullfile(here, '..'); res = @(varargin) fullfile(root, 'results', varargin{:});

fprintf('\n=== Exp. 10: nonlinear memory and chance baselines (Exp. 1 runs) ===\n');
old = cd(root); T = collect_grid('exp10_nonlinear'); cd(old);
g = @(f) arrayfun(@(t) double(t.(f)), T)';
pt = reshape([T.pt], 3, [])'; MC = g('MC'); NLC = g('NLC'); MCs = g('MC_shuf'); NLs = g('NLC_shuf');
lc = g('lam_cond'); q = g('rateE') < 4; sd = pt(:,3);
E1 = load_scan('exp1_regime');
dd = nan(size(MC)); st = false(size(MC));
for i = 1:numel(MC), j = find(all(abs(E1.pt - pt(i,:)) < 1e-9, 2)); dd(i) = MC(i) - E1.MC(j); st(i) = (E1.rateE(j) < 4) ~= q(i); end
k = abs(dd) > 1e-9;
fprintf('repeat of the Exp. 1 driven runs: identical %d of %d; differing %d (lam_cond>0 in %d), |dMC| median %.3f max %.3f; state changed %d\n', ...
    sum(~k), numel(k), sum(k), sum(k & E1.lam_cond(arrayfun(@(i) find(all(abs(E1.pt - pt(i,:)) < 1e-9, 2)), 1:numel(MC))) > 0), median(abs(dd(k))), max(abs(dd(k))), sum(st));
fprintf('median quiet-state rate in Exp. 1: %.2f Hz\n', median(E1.rateE(E1.rateE < 4)));
fprintf('chance MC %.3f +- %.3f (max %.3f); chance NLC %.3f +- %.3f (max %.3f)\n', mean(MCs), std(MCs), max(MCs), mean(NLs), std(NLs), max(NLs));
fprintf('NLC quiet %.2f +- %.2f, active %.2f +- %.2f\n', mean(NLC(q)), std(NLC(q)), mean(NLC(~q)), std(NLC(~q)));
[d, ci] = cboot(@(i) mean(NLC(i(q(i)))) - mean(NLC(i(~q(i)))), sd); fprintf('NLC quiet - active %.2f [%.2f, %.2f]\n', d, ci);
[r, ci] = cboot(@(i) corr_rank(NLC(i), lc(i)), sd); fprintf('Spearman NLC vs lam_cond %.2f [%.2f, %.2f]\n', r, ci);
[r, ci] = cboot(@(i) corr_rank(NLC(i), MC(i)), sd); fprintf('Spearman NLC vs MC %.2f [%.2f, %.2f]\n', r, ci);
fprintf('NLC > MC in %.0f%% of quiet and %.0f%% of active runs\n', 100*mean(NLC(q) > MC(q)), 100*mean(NLC(~q) > MC(~q)));
[m, im] = max(NLC); fprintf('max NLC %.2f at lam_cond %.2f (quiet %d)\n', m, lc(im), q(im));
ed = [-Inf -10 -5 -2 0 Inf];
for b = 1:numel(ed)-1, s = lc >= ed(b) & lc < ed(b+1); fprintf('  lam_cond [%g,%g): NLC %.2f (n=%d)\n', ed(b), ed(b+1), mean(NLC(s)), sum(s)); end
K = size(E1.MCk, 2); qa = E1.rateE < 4;
mq = mean(E1.MCk(qa,:)); ma = mean(E1.MCk(~qa,:));
fprintf('MC_k quiet lags 1-5: %s; active: %s\n', mat2str(mq(1:5), 2), mat2str(ma(1:5), 2));
fprintf('lag where mean MC_k < 0.1: quiet %d ms, active %d ms\n', 20*find(mq < 0.1, 1), 20*find(ma < 0.1, 1));
fprintf('share of MC at lags > 100 ms: quiet %.2f, active %.2f\n', sum(mq(6:end))/sum(mq), sum(ma(6:end))/sum(ma));

fprintf('\n=== Memory versus distance to own linear edge (connectivities 1-3) ===\n');
[dd, mc, qq, grp] = distance_data(root);
s = dd > 0; fprintf('runs below the edge %d (quiet %d, active %d); above %d (quiet %d)\n', sum(s), sum(s & qq), sum(s & ~qq), sum(~s), sum(~s & qq));
fprintf('Spearman MC vs d, all runs %.2f (n=%d); quiet runs below edge %.2f (n=%d); active runs below edge %.2f (n=%d)\n', ...
    corr_rank(mc, dd), numel(mc), corr_rank(mc(s & qq), dd(s & qq)), sum(s & qq), corr_rank(mc(s & ~qq), dd(s & ~qq)), sum(s & ~qq));
for b = [0 0.05; 0.05 0.15; 0.15 0.3; 0.3 1]'
    t = qq & dd >= b(1) & dd < b(2); fprintf('  quiet, d in [%.2f,%.2f): MC %.2f (n=%d)\n', b, mean(mc(t)), sum(t));
end
t = ~qq & dd > 0; fprintf('  active below the edge: MC %.2f (n=%d)\n', mean(mc(t)), sum(t));
% rate-matched comparison: quiet runs within a narrow rate band
r1 = [E1.rateE]; t1 = qa & r1 > 1.2 & r1 < 2.5;
fprintf('quiet runs with rates 1.2-2.5 Hz: n=%d, Spearman MC vs rate %.2f, MC vs lam_cond %.2f\n', sum(t1), corr_rank(E1.MC(t1), r1(t1)), corr_rank(E1.MC(t1), E1.lam_cond(t1)));

fprintf('\n=== Exp. 9: numerical checks ===\n');
L = load(res('exp9_numerics', 'numerics.mat'));
for r = 1:numel(L.res)
    x = L.res(r);
    fprintf('%-14s ref: rate %5.2f MC %4.2f lamA %6.2f lamC %6.2f | 24-s lamA %6.2f +- %4.2f | dt/2: rate %5.2f MC %4.2f lamA %6.2f lamC %6.2f | long sets: rate %5.2f MC %4.2f chance %.3f\n', ...
        mat2str(x.pt, 3), x.ref.rateE, x.ref.MC, x.ref.lam_aut, x.ref.lam_cond, mean(x.lam_blocks), std(x.lam_blocks), ...
        x.dt.rateE, x.dt.MC, x.dt.lam_aut, x.dt.lam_cond, x.mem.rateE, x.mem.MC, x.mem.MC_shuf);
end
same = arrayfun(@(x) (x.ref.rateE < 4) == (x.dt.rateE < 4), L.res);
far = arrayfun(@(x) abs(x.ref.lam_aut) > 5, L.res);
fprintf('state unchanged with dt/2 in %d of %d points\n', sum(same), numel(same));
dl = arrayfun(@(x) abs(x.dt.lam_cond - x.ref.lam_cond), L.res(same));
dm = arrayfun(@(x) abs(x.dt.MC - x.ref.MC), L.res(same));
fprintf('same-state points: |d lam_cond| median %.2f max %.2f; |d MC| median %.2f max %.2f\n', median(dl), max(dl), median(dm), max(dm));
dla = arrayfun(@(x) mean(x.lam_blocks) - x.ref.lam_aut, L.res);
fprintf('24-s minus 4-s lam_aut: points with |lam_aut|>5: %s; near zero: %s\n', mat2str(dla(far), 2), mat2str(dla(~far), 2));
fprintf('block s.d. of lam_aut: median %.2f, max %.2f\n', median(arrayfun(@(x) std(x.lam_blocks), L.res)), max(arrayfun(@(x) std(x.lam_blocks), L.res)));
rl = arrayfun(@(x) x.mem.MC / x.ref.MC, L.res); sm = arrayfun(@(x) (x.ref.rateE < 4) == (x.mem.rateE < 4), L.res);
fprintf('long-set MC / standard MC (same state): median %.2f, range %.2f-%.2f; chance %.3f-%.3f\n', median(rl(sm)), min(rl(sm)), max(rl(sm)), min(arrayfun(@(x) x.mem.MC_shuf, L.res)), max(arrayfun(@(x) x.mem.MC_shuf, L.res)));
for k = 1:3
    b = L.base(k); fprintf('no recurrence seed %d: same drive MC %.3f NLC %.3f (%.2f Hz); matched MC %.3f NLC %.3f (%.2f Hz)\n', k, b.same_drive.MC, b.same_drive.NLC, b.same_drive.rateE, b.matched.MC, b.matched.NLC, b.matched.rateE);
end
f = res('exp9_numerics', 'replicates.mat');
if exist(f, 'file')
    R = load(f); R = R.rep;
    for r = 1:numel(R)
        fprintf('%-14s 4 realisations: active fraction dt 0.5 %.2f / 0.25 %.2f; MC %.2f / %.2f; lamC %.2f / %.2f\n', mat2str(R(r).pt, 3), ...
            mean(R(r).rate(1,:) >= 4), mean(R(r).rate(2,:) >= 4), mean(R(r).MC(1,:)), mean(R(r).MC(2,:)), mean(R(r).lam_cond(1,:)), mean(R(r).lam_cond(2,:)));
    end
    af = arrayfun(@(x) abs(mean(x.rate(1,:) >= 4) - mean(x.rate(2,:) >= 4)), R);
    fprintf('|difference in active fraction| between steps: median %.2f, max %.2f\n', median(af), max(af));
    mix = arrayfun(@(x) any(x.rate(1,:) >= 4) && any(x.rate(1,:) < 4), R);
    fprintf('points whose state depends on the realisation (dt 0.5): %d of %d\n', sum(mix), numel(R));
end

% long run at the point with the largest long/standard discrepancy: MC in standard-length windows
old = cd(root); addpath(fullfile(root, 'src')); p = [2.6 0.25 3]; P = set_composition(base_params(), p(2), p(1), -55); C = build_connectivity(P, 3);
randn('state', 1003); u = randn(7400, 1);
od = simulate_net(P, C, struct('u', u, 'noise', true, 'lyap', true, 'lyap_skip', 200, 'seed', 2003)); cd(old);
w = arrayfun(@(w0) memory_capacity(od.R(:, w0 + (1:2000)), u(w0 + (1:2000)), 200, 1200, 50), [0 1800 3600 5400]);
fprintf('point %s: MC of the long run over the long split %.2f; in four standard windows %s\n', mat2str(p, 3), memory_capacity(od.R, u, 200, 4800, 50), mat2str(w, 2));
lq = arrayfun(@(x) x.mem.MC, L.res); sq = arrayfun(@(x) x.mem.rateE < 4, L.res);
fprintf('long-set MC: quiet %.2f (n=%d), active %.2f (n=%d)\n', mean(lq(sq)), sum(sq), mean(lq(~sq)), sum(~sq));

fprintf('\n=== Exp. 3c: arclength continuation ===\n');
for s = 1:3, for fN = [0 0.25 0.5 0.75]
    A = load(res('exp3c_arclength', sprintf('seed%d_fN%02d.mat', s, round(100*fN))));
    B = load(res('exp3b_continuation', sprintf('seed%d_fN%02d_Mg100.mat', s, round(100*fN))));
    o = A.out; i1 = find(diff(o(:,1)) < 0, 1);
    fprintf('seed %d fN %.2f: fold %.4f (Newton branch end %.4f), unstable eigenvalues just after the fold %d, first loss %.4f (Newton %.4f), last stable low-rate %.4f\n', ...
        s, fN, A.g_fold, B.gE_end, o(min(i1+1, end), 6), A.g_loss, pick(B), A.g_lastq);
end, end

A = load(res('exp3c_arclength', 'seed1_fN00.mat')); o = A.out; stb = o(:,6) == 0 & o(:,3) < 4;
d = diff([0; stb; 0]); a1 = find(d == 1); b1 = find(d == -1) - 1;
for k = 1:numel(a1), fprintf('  seed 1 fN 0 stable low-rate segment %d: gE %.4f-%.4f, rate %.2f-%.2f Hz\n', k, min(o(a1(k):b1(k),1)), max(o(a1(k):b1(k),1)), min(o(a1(k):b1(k),3)), max(o(a1(k):b1(k),3))); end
% bistable points: distance below the edge
C0 = load(res('exp12_coexistence', 'coexistence.mat')); bi0 = C0.out(:,1) < 4 & C0.out(:,2) < 1e-2 & C0.out(:,3) >= 4;
fprintf('bistable points lie %s %% below the edge\n', mat2str(round(100 * (1 - C0.rows(bi0,1) ./ C0.rows(bi0,4)))', 2));

fprintf('\n=== Exp. 12: coexistence below the linear edge ===\n');
C = load(res('exp12_coexistence', 'coexistence.mat')); o = C.out;
bi = o(:,1) < 4 & o(:,2) < 1e-2 & o(:,3) >= 4;   % quiet fixed point (converging) and active attractor
fprintf('points %d; near quiet FP: returned to it %d; from driven active state: active %d; bistable %d; default start active %d\n', ...
    size(C.rows, 1), sum(o(:,1) < 4 & o(:,2) < 1e-2), sum(o(:,3) >= 4), sum(bi), sum(o(:,5) >= 4));
fprintf('active attractor: rates %.1f-%.1f Hz, fixed point in %d\n', min(o(o(:,3)>=4,3)), max(o(o(:,3)>=4,3)), sum(o(:,3) >= 4 & o(:,4) < 1e-4));
for f = [0 0.25 0.5 0.75], s = C.rows(:,2) == f; fprintf('  fN %.2f: %d points, bistable %d\n', f, sum(s), sum(bi & s)); end

fprintf('\n=== Exp. 13: slowest mode ===\n');
S = load(res('exp13_slow_mode', 'slow_mode.mat')); o = S.out;
% columns: seed gE tauN decay imp imp20 imp60 imp100 imp200 dvar iso resid
fprintf('n=%d fixed points; Lyapunov residual max %.1e\n', size(o,1), max(o(:,12)));
fprintf('slowest decay x tau_N: %s\n', mat2str((o(:,4) .* o(:,3) / 1000)', 2));
fprintf('slow-mode share of impulse-response power (20-1000 ms): median %.3f, range %.3f-%.3f\n', median(o(:,5)), min(o(:,5)), max(o(:,5)));
for t = [25 50 100 200], k = o(:,3) == t; fprintf('  tau_N %3d: %s\n', t, mat2str(o(k,5)', 2)); end
for sg = unique(o(:,1:2), 'rows')', k = o(:,1) == sg(1) & o(:,2) == sg(2); fprintf('  seed %d gE %.1f: share vs tau_N %s, monotonic decrease %d\n', sg(1), sg(2), mat2str(o(k,5)', 2), all(diff(o(k,5)) < 0)); end
fprintf('single-lag ratio at 200 ms: %s (can exceed 1: interference between modes)\n', mat2str(o(:,9)', 2));
fprintf('change of stationary rate variance without the slow mode: %.3f to %.3f; isolated-mode measure median %.1e\n', min(o(:,10)), max(o(:,10)), median(o(:,11)));
E1b = load_scan('exp1_regime'); qb = E1b.rateE < 4; mk = mean(E1b.MCk(qb,:));
fprintf('mean quiet MC_k at lags 100 and 200 ms: %.3f, %.3f\n', mk(5), mk(10));

fprintf('\n=== Exp. 8: homeostatic rules ===\n');
v = {'eta0.3', 'eta0.1', 'eta0.03', 'tanh', 'cont', 'floor'};
for k = 1:numel(v)
    for c = 1:2, for s = 1:3
        fn = res('exp8_homeostasis', sprintf('seed%d_cond%d_%s.mat', s, c, v{k}));
        if ~exist(fn, 'file'), continue; end
        L = load(fn); H = L.H; n = numel(H.rateE); i0 = round(n/3) + 1; tsec = L.h.T_epoch * 0.02;
        act = H.rateE > 4; ep = sum(diff([0; act(:)]) == 1);
        on = find(diff([0; act(:)]) == 1); off = find(diff([act(:); 0]) == -1); dur = (off - on + 1) * tsec;
        quiet = ~act; qr = H.rateE(quiet); tmin = 0.25; if isfield(L.h, 'theta_min'), tmin = L.h.theta_min; end
        fprintf('%-7s c%d s%d: %4d epochs (%.0f s each) | episodes %3d, duration %s s (median %.0f), active fraction (last 2/3) %.2f | mean rate (last 2/3) %.2f Hz = %.2f x target | quiet epochs: median rate %.2f x target, max %.2f Hz (target %.2f) | theta %.2f-%.2f, at bound %d\n', ...
            v{k}, c, s, n, tsec, ep, rng(dur), median([dur; NaN](1:max(1,numel(dur)))), mean(act(i0:end)), mean(H.rateE(i0:end)), mean(H.rateE(i0:end)) / L.target, ...
            median(qr) / L.target, max(qr), L.target, min(H.theta), max(H.theta), sum(H.theta <= tmin + 1e-9 | H.theta >= 6 - 1e-9));
    end, end
end

fprintf('\n=== Exp. 11: restoration from actual states ===\n');
for c = 1:2, for s = 1:3
    fn = res('exp11_restore', sprintf('seed%d_cond%d.mat', s, c));
    if ~exist(fn, 'file'), continue; end
    L = load(fn); nm = fieldnames(L.R);
    O = load(res('exp6_vision', sprintf('seed%d_cond%d.mat', s, c)));
    fprintf('c%d s%d: reproduces Exp. 6 theta trajectory: max |d theta| %.1e\n', c, s, max(abs(L.H.theta - O.Out.hist_dep.theta)));
    for q = 1:numel(nm)
        X = L.R.(nm{q}); th = [X.theta]; rt = [X.rateE]; mcv = [X.MC]; lcv = [X.lam_cond];
        fprintf('   %-10s theta %s\n              rate  %s\n              MC    %s\n              hyperactive (>=4 Hz) %d/%d; MC %.2f-%.2f\n', nm{q}, mat2str(th, 2), mat2str(rt, 2), mat2str(mcv, 2), sum(rt >= 4), numel(rt), min(mcv), max(mcv));
        qth = th(rt < 4); hth = th(rt >= 4);
        fprintf('              largest theta with quiet outcome %.2f, smallest theta with hyperactive outcome %.2f; MC quiet outcomes %.2f-%.2f, hyperactive %.2f-%.2f\n', ...
            max([qth NaN]), min([hth NaN]), min([mcv(rt < 4) NaN]), max([mcv(rt < 4) NaN]), min([mcv(rt >= 4) NaN]), max([mcv(rt >= 4) NaN]));
        if q == 1, fprintf('              lam_cond of hyperactive outcomes %s\n', mat2str(lcv(rt >= 4), 2)); end
        if q > 1
            F = L.R.full; rf = [F.rateE];
            fprintf('              rescued (full hyperactive -> quiet) %d of %d hyperactive epochs; rate reduction in still-hyperactive epochs %s\n', ...
                sum(rf >= 4 & rt < 4), sum(rf >= 4), mat2str(rt(rf >= 4 & rt >= 4) ./ rf(rf >= 4 & rt >= 4), 2));
        end
    end
end, end
end

function s = rng(x)
if isempty(x), s = '-'; else, s = sprintf('%.0f-%.0f', min(x), max(x)); end
end

function e = pick(B)
if isfinite(B.gE_unstable(3)), e = B.gE_unstable(3); else, e = B.gE_end; end
end

function [d, mc, q, grp] = distance_data(root)
rd = fullfile(root, 'results', 'exp3b_continuation');
ld = @(s, f, m) load(fullfile(rd, sprintf('seed%d_fN%02d_Mg%03d.mat', s, round(100*f), round(100*m))));
ed = @(L, t) edge_t(L, t);
d = []; mc = []; q = []; grp = [];
E1 = load_scan('exp1_regime');
for s = 1:3, for f = [0 0.25 0.5 0.75]
    e = ed(ld(s, f, 1), 3); i = E1.pt(:,3) == s & E1.pt(:,2) == f;
    d = [d; 1 - E1.pt(i,1) / e]; mc = [mc; E1.MC(i)]; q = [q; E1.rateE(i) < 4]; grp = [grp; ones(sum(i), 1)]; %#ok<AGROW>
end, end
E2 = load_scan('exp2_kinetics'); tN = [25 50 100 200];
for s = 1:2, for t = [1 2 4]
    e = ed(ld(s, 0.5, 1), t); i = E2.pt(:,3) == s & E2.pt(:,2) == tN(t);
    d = [d; 1 - E2.pt(i,1) / e]; mc = [mc; E2.MC(i)]; q = [q; E2.rateE(i) < 4]; grp = [grp; 2 * ones(sum(i), 1)]; %#ok<AGROW>
end, end
E3 = load_scan('exp3_magnesium');
for s = 1:2, for m = [0.5 2]
    e = ed(ld(s, 0.5, m), 3); i = E3.pt(:,3) == s & abs(E3.pt(:,2) - m) < 1e-9;
    d = [d; 1 - E3.pt(i,1) / e]; mc = [mc; E3.MC(i)]; q = [q; E3.rateE(i) < 4]; grp = [grp; 3 * ones(sum(i), 1)]; %#ok<AGROW>
end, end
q = logical(q);
end

function e = edge_t(L, t)
if isfinite(L.gE_unstable(t)), e = L.gE_unstable(t);
elseif isfield(L, 'gE_end') && isfinite(L.gE_end), e = L.gE_end;
else, e = L.gEs(end); end
end

function [est, ci] = cboot(fun, seed, nb)
if nargin < 3, nb = 2000; end
seed = seed(:); us = unique(seed); est = fun((1:numel(seed))');
st = rand('state'); rand('state', 2); v = nan(nb, 1);
for b = 1:nb
    pick = us(ceil(rand(numel(us), 1) * numel(us)));
    idx = cell2mat(arrayfun(@(s) find(seed == s), pick, 'UniformOutput', false));
    v(b) = fun(idx);
end
rand('state', st); v = sort(v(isfinite(v)));
ci = [v(max(1, round(0.025*numel(v)))), v(round(0.975*numel(v)))];
end

function r = corr_rank(a, b)
a = a(:); b = b(:); ok = isfinite(a) & isfinite(b); a = a(ok); b = b(ok);
[~, ia] = sort(a); ra(ia) = 1:numel(a); [~, ib] = sort(b); rb(ib) = 1:numel(b);
c = corrcoef(ra, rb); r = c(1,2);
end
