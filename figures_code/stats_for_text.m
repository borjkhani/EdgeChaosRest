function stats_for_text()
%STATS_FOR_TEXT  Print every number quoted in the manuscript text.
%   Pooled statistics use the distinct runs of Experiments 1-3: the
%   reference conditions of Experiment 2 (tau_N = 100 ms) and Experiment 3
%   ([Mg2+] = 1 mM) are identical to Experiment 1 at fN = 0.5, seeds 1-2,
%   and are therefore counted once. Confidence intervals resample
%   connectivities (cluster bootstrap), because runs that share a
%   connectivity are not independent.
here = fileparts(mfilename('fullpath')); addpath(here);
addpath(fullfile(here, '..', 'src'));
root = fullfile(here, '..');
E1 = load_scan('exp1_regime'); E2 = load_scan('exp2_kinetics'); E3 = load_scan('exp3_magnesium');
k2 = E2.pt(:,2) ~= 100; k3 = abs(E3.pt(:,2) - 1) > 1e-9;
f = @(n) [E1.(n); E2.(n)(k2); E3.(n)(k3)];
rate = f('rateE'); MC = f('MC'); lc = f('lam_cond'); la = f('lam_aut');
rk = f('rhoK'); fp = f('has_fp'); q = rate < 4; sdp = f2(E1, E2, E3, k2, k3);
fprintf('\n=== Pooled, Experiments 1-3 (distinct runs) ===\n');
fprintf('runs %d (quiet %d, active %d)\n', numel(rate), sum(q), sum(~q));
fprintf('quiet rates %.2f-%.2f Hz; active rates %.2f-%.2f Hz\n', min(rate(q)), max(rate(q)), min(rate(~q)), max(rate(~q)));
fprintf('fixed points %d (quiet %d, active %d); max rhoK quiet %.3f, overall %.3f; all < 1: %d\n', ...
    sum(fp==1), sum(fp==1 & q), sum(fp==1 & ~q), max(rk(fp==1 & q)), max(rk(fp==1)), all(rk(fp==1) < 1));
fprintf('MC quiet %.2f +- %.2f (max %.2f) | active %.2f +- %.2f\n', mean(MC(q)), std(MC(q)), max(MC(q)), mean(MC(~q)), std(MC(~q)));
edges = [-Inf -20 -10 -5 -2 0 Inf];
for b = 1:numel(edges)-1
    s = lc >= edges(b) & lc < edges(b+1);
    fprintf('  lam_cond in [%g,%g): MC %.2f (n=%d)\n', edges(b), edges(b+1), mean(MC(s)), sum(s));
end
[m, im] = max(MC); fprintf('max MC %.2f at lam_cond %.2f (rate %.2f Hz)\n', m, lc(im), rate(im));
[d, dci] = cboot(@(i) mean(MC(i(q(i)))) - mean(MC(i(~q(i)))), sdp);
fprintf('MC quiet - active: %.2f, 95%% CI [%.2f, %.2f] (cluster bootstrap over %d connectivities)\n', d, dci, numel(unique(sdp)));
[r, rci] = cboot(@(i) corr_rank(MC(i), lc(i)), sdp);
fprintf('Spearman MC vs lam_cond, all runs: %.2f, 95%% CI [%.2f, %.2f]\n', r, rci);
[r, rci] = cboot(@(i) corr_rank(MC(i(q(i))), lc(i(q(i)))), sdp);
fprintf('Spearman MC vs lam_cond, quiet runs: %.2f, 95%% CI [%.2f, %.2f] (n=%d)\n', r, rci, sum(q));
sel = fp==1 & q; fprintf('Spearman MC vs rhoK, quiet fixed points (n=%d): %.2f\n', sum(sel), corr_rank(MC(sel), rk(sel)));
% stable fixed points with rhoK > 1: leading eigenvalue of K
for E = {E1}, X = E{1}; ii = find(X.has_fp == 1 & X.rhoK > 1);
    for i = ii', Pp = set_composition(base_params(), X.pt(i,2), X.pt(i,1), -55); Cc = build_connectivity(Pp, X.pt(i,3));
        x = relax_to_fp(Pp, Cc, [], 250); [Jj, inf] = jacobian_fp(Pp, Cc, x); Bb = mg_block(x.V, Pp);
        rowE = Pp.gA * inf.dfE + Pp.gN * Bb .* inf.dfN; Dd = diag(inf.dnu);
        Kk = [diag(rowE ./ inf.Gloc) * Cc.JE * Dd(Cc.isE, Cc.isE), diag(Pp.gG * inf.dfI ./ inf.Gloc) * Cc.JI * Dd(~Cc.isE, ~Cc.isE)];
        ek = eig(Kk); [~, m1] = max(real(ek)); er = real(ek(abs(imag(ek)) < 1e-9));
        fprintf('rhoK>1 fixed point %s (rate %.1f Hz): lead eig of K %.3f%+.3fi, largest real eig %.2f, max Re Jac %.2f /s\n', mat2str(X.pt(i,:), 3), X.rateE(i), real(ek(m1)), imag(ek(m1)), max(er), 1000 * max(real(eig(Jj))));
    end
end
ra = []; rd1 = []; hf = []; old1 = cd(root); nms = {'exp1_regime', 'exp2_kinetics', 'exp3_magnesium'};
for k = 1:3
    T1 = collect_grid(nms{k});
    for t = T1
        if (k == 2 && t.pt(2) == 100) || (k == 3 && abs(t.pt(2) - 1) < 1e-9), continue; end
        ra(end+1) = double(t.rateE_aut); rd1(end+1) = double(t.rateE); hf(end+1) = double(t.has_fp); %#ok<AGROW>
    end
end
cd(old1); s1 = hf == 1 & rd1 >= 4;
fprintf('active runs on fixed points (distinct): %d; fixed-point rate >= 4 Hz: %d (%.1f-%.1f Hz); < 4 Hz: %d\n', sum(s1), sum(s1 & ra >= 4), min(ra(s1 & ra >= 4)), max(ra(s1 & ra >= 4)), sum(s1 & ra < 4));
gs = unique(E1.pt(:,1))'; for fN = [0 0.25 0.5 0.75], mm = arrayfun(@(g) mean(E1.MC(E1.pt(:,1) == g & E1.pt(:,2) == fN)), gs); fprintf('fN %.2f mean MC (all runs) by gE: %s\n', fN, mat2str(mm, 2)); end
act = ~q; v = sort(la(act)); fprintf('active lam_aut 10-90%%: %.1f..%.1f; active runs on stable fixed points: %d\n', v(round(0.1*end)), v(round(0.9*end)), sum(fp==1 & act));
fprintf('lam_cond - lam_aut: quiet median %.2f; active median %.2f, frac positive %.2f (n=%d)\n', ...
    median(lc(q)-la(q)), median(lc(act)-la(act)), mean(lc(act) > la(act)), sum(act));

fprintf('\n=== Experiment 1 (%d runs, %d connectivities) ===\n', numel(E1.MC), numel(unique(E1.pt(:,3))));
a1 = E1.rateE >= 4; gs = unique(E1.pt(:,1))';
for fN = [0 0.25 0.5 0.75]
    s = E1.pt(:,2)==fN & E1.rateE<4; sa = E1.pt(:,2)==fN;
    [g50, ci] = logistic_edge(E1.pt(sa,1), a1(sa), E1.pt(sa,3));
    na = arrayfun(@(g) sum(a1(sa & E1.pt(:,1)==g)), gs);
    fprintf('fN=%.2f: largest quiet gE %.1f; quiet MC at gE=1.2 %.2f; P(active)=0.5 at gE %.2f [%.2f, %.2f]; active per gE %s\n', ...
        fN, max(E1.pt(s,1)), mean(E1.MC(s & E1.pt(:,1)==1.2)), g50, ci, mat2str(na));
    fprintf('        quiet runs %d/%d, quiet MC %.2f, active MC %.2f, all-run MC %.2f\n', sum(s), sum(sa), mean(E1.MC(s)), mean(E1.MC(sa & a1)), mean(E1.MC(sa)));
end
% shift of the 50% point between fN = 0 and fN > 0, paired over connectivities
for fN = [0.25 0.5 0.75]
    [d, dci] = cboot(@(i) g50_of(E1, a1, i, fN) - g50_of(E1, a1, i, 0), E1.pt(:,3));
    fprintf('shift of P(active)=0.5 point, fN=%.2f vs 0: %.2f [%.2f, %.2f]\n', fN, d, dci);
end
nq = arrayfun(@(sd) all(~a1(E1.pt(:,3)==sd & E1.pt(:,2)==0)), unique(E1.pt(:,3)));
fprintf('connectivities quiet at all couplings (fN=0): %d\n', sum(nq));
for sd = unique(E1.pt(:,3))'
    s = E1.pt(:,3)==sd; fprintf('  seed %2d: first active gE by fN %s\n', sd, mat2str(arrayfun(@(fN) min([E1.pt(s & E1.pt(:,2)==fN & a1, 1); Inf]), [0 0.25 0.5 0.75]), 2));
end
s0 = E1.pt(:,2)==0 & E1.rateE<4; s1 = E1.pt(:,2)>0 & E1.rateE<4; v = sort(E1.lam_aut(s1));
fprintf('quiet lam_aut fN=0: [%.1f, %.1f]; fN>0: 10-90%% [%.1f, %.1f], full [%.1f, %.1f], n=%d, n>=0: %d\n', ...
    min(E1.lam_aut(s0)), max(E1.lam_aut(s0)), v(round(0.1*end)), v(round(0.9*end)), min(v), max(v), numel(v), sum(v>=0));
r = []; nr = [];
for sd = unique(E1.pt(:,3))', for fN = [0 0.25 0.5 0.75]
    s = E1.pt(:,3)==sd & E1.pt(:,2)==fN & E1.rateE<4;
    if sum(s) >= 3, c = corrcoef(E1.pt(s,1), E1.MC(s)); r(end+1) = c(1,2); nr(end+1) = sum(s); %#ok<AGROW>
        fprintf('  seed %d fN %.2f: n=%d, r(gE,MC)=%.2f, MC %.2f-%.2f\n', sd, fN, sum(s), c(1,2), min(E1.MC(s)), max(E1.MC(s))); end
end, end
fprintf('quiet MC vs gE within connectivity: %d groups, r>0 in %d, median r %.2f; groups with >=6 quiet couplings: %d, r>0 in %d\n', ...
    numel(r), sum(r>0), median(r), sum(nr>=6), sum(r(nr>=6)>0));

fprintf('\n=== Experiment 2 ===\n');
d = 0; n = 0;
for g = unique(E2.pt(:,1))', for sd = 1:2
    s = E2.pt(:,1)==g & E2.pt(:,3)==sd & E2.has_fp==1;
    if sum(s) > 1, d = max(d, max(E2.rhoK(s)) - min(E2.rhoK(s))); n = n + 1; end
end, end
fprintf('max spread of rhoK across tau_N: %.2e over %d (gE, seed) groups\n', d, n);
for t = [25 50 100 200]
    s = E2.pt(:,2)==t & E2.rateE<4;
    fprintf('tau_N=%3d: quiet %d, lam_aut %.2f (x tau_N %.2f), lam_cond %.2f, MC %.2f, largest quiet gE %.1f\n', ...
        t, sum(s), mean(E2.lam_aut(s)), mean(E2.lam_aut(s))*t/1000, mean(E2.lam_cond(s)), mean(E2.MC(s)), max(E2.pt(s,1)));
end

fprintf('\n=== Experiment 3 ===\n');
for m = [0.25 0.5 1 2]
    s = abs(E3.pt(:,2)-m)<1e-9; sq = s & E3.rateE<4; sf = s & E3.has_fp==1;
    gq = E3.pt(sq,1); if isempty(gq), gq = NaN; end
    fprintf('Mg=%.2f: quiet %d/%d, largest quiet gE %.1f, quiet MC %.2f, active rate %.1f Hz, NSC dRhoK median %.3f max %.3f\n', ...
        m, sum(sq), sum(s), max(gq), mean(E3.MC(sq)), mean(E3.rateE(s & ~sq)), median(E3.rhoK(sf)-E3.rhoK_noNSC(sf)), max(E3.rhoK(sf)-E3.rhoK_noNSC(sf)));
end

dn = [];
for E = {E1, E2, E3}, X = E{1}; sN = X.has_fp == 1 & ~(X.pt(:,2) == 0 & numel(unique(X.pt(:,2))) == 4); dn = [dn; X.rhoK(sN) - X.rhoK_noNSC(sN)]; end %#ok<AGROW>
fprintf('NSC effect on rhoK at NMDA fixed points (Exp. 1-3, duplicates included): min %.4f, median %.3f, n<=0: %d\n', min(dn), median(dn), sum(dn <= 0));

fprintf('\n=== Experiment 4 (heterogeneous tau_N; pt = [gE fN het seed]) ===\n');
E4 = load_scan('exp4_tau_heterogeneity'); dMC = []; rl = []; sw = 0;
for sd = 1:3, for fN = [0.5 0.75], for g = [1.2 1.4]
    i0 = find(E4.pt(:,1)==g & E4.pt(:,2)==fN & E4.pt(:,3)==0 & E4.pt(:,4)==sd);
    i1 = find(E4.pt(:,1)==g & E4.pt(:,2)==fN & E4.pt(:,3)==1 & E4.pt(:,4)==sd);
    if E4.rateE(i0) < 4 && E4.rateE(i1) < 4, dMC(end+1) = E4.MC(i1)-E4.MC(i0); rl(end+1) = E4.lam_cond(i1)/E4.lam_cond(i0); %#ok<AGROW>
    elseif E4.rateE(i0) < 4, sw = sw + 1; fprintf('  switched to active: gE %.1f fN %.2f seed %d (%.2f -> %.2f Hz)\n', g, fN, sd, E4.rateE(i0), E4.rateE(i1)); end
end, end, end
fprintf('quiet pairs %d: dMC mean %.3f range [%.3f, %.3f]; lam_cond ratio %.2f; switched %d\n', numel(dMC), mean(dMC), min(dMC), max(dMC), mean(rl), sw);
rand('state', 100); tt = exp(log(30) + rand(1e6,1)*(log(330)-log(30))); fprintf('log-uniform [30,330] ms: mean %.0f ms\n', mean(tt));

fprintf('\n=== Experiment 5 (hysteresis) ===\n');
for c = 1:4, for sd = 1:2
    L = load(fullfile(root, 'results', 'exp5_hysteresis', sprintf('cond%d_seed%d.mat', c, sd)));
    i = find(L.rate(:,1,1) > 4, 1); j = find(L.rate(:,2,2) > 4, 1); jl = find(L.rate(:,2,1) > 4, 1);
    if isempty(jl), mlow = max(L.rate(:,2,1)); jls = 'none'; else, mlow = max(L.rate(1:jl-1,2,1)); jls = sprintf('%.1f', L.theta(jl)); end
    bis = sum(abs(L.rate(:,2,1) - L.rate(:,2,2)) > 2);
    fprintf('cond%d s%d: intact jump %.1f-%.1f; degenerated: high branch from %.1f, low branch jumps at %s, low-branch max %.2f Hz, bistable grid points %d\n', ...
        c, sd, L.theta(i-1), L.theta(i), L.theta(j), jls, mlow, bis);
end, end

fprintf('\n=== Experiment 6 (vision protocol) ===\n');
for c = 1:2, for sd = 1:3
    L = load(fullfile(root, 'results', 'exp6_vision', sprintf('seed%d_cond%d.mat', sd, c))); O = L.Out; H = O.hist_dep;
    ep = H.rateE > 4; tgt = mean(O.r_target(1:160)); R2 = O.rest(2); R1 = O.rest(1);
    kfull = find(abs(R2.hist_rec.rateE - tgt) < 0.2*tgt, 1);
    j1 = find(H.rateE > 4, 1);
    fprintf('cond%d s%d: target %.2f Hz | episodes %d (%.1f-%.1f Hz), theta %.2f-%.2f, theta at first episode %.2f, theta_dep %.2f\n', ...
        c, sd, tgt, sum(ep), min(H.rateE(ep)), max(H.rateE(ep)), min(H.theta), max(H.theta), H.theta(max(j1-1,1)), O.theta_dep);
    fprintf('   healthy MC %.2f | full restoration: rate %.2f MC %.2f lamC %.2f | +Mg MC %.2f rate %.2f | +blk MC %.2f rate %.2f | rescaled MC %.2f (recovery epoch %d)\n', ...
        O.healthy.MC, R2.frozen.rateE, R2.frozen.MC, R2.frozen.lam_cond, R2.high_mg.MC, R2.high_mg.rateE, R2.nmda_block.MC, R2.nmda_block.rateE, R2.rescaled.MC, kfull);
    fprintf('   partial restoration: rate %.2f MC %.2f | re-scaling under partial input: episodes %d\n', R1.frozen.rateE, R1.frozen.MC, sum(R1.hist_rec.rateE > 4));
end, end

fprintf('\n=== Experiment 7 (current-based AMPA/GABA) ===\n');
E7 = load_scan('exp7_current_based');
[m, im] = max(E7.MC);
fprintf('rates %.1f-%.1f Hz; lam_cond<lam_aut in %.0f%%; max MC %.2f at lam_cond %.2f (rate %.1f Hz)\n', ...
    min(E7.rateE), max(E7.rateE), 100*mean(E7.lam_cond < E7.lam_aut), m, E7.lam_cond(im), E7.rateE(im));
fprintf('MC lam_cond>0 %.2f (n=%d); lam_cond<0 %.2f (n=%d); fixed points %d, max rhoK %.3f\n', ...
    mean(E7.MC(E7.lam_cond>0)), sum(E7.lam_cond>0), mean(E7.MC(E7.lam_cond<0)), sum(E7.lam_cond<0), sum(E7.has_fp==1), max(E7.rhoK(E7.has_fp==1)));

fprintf('\n=== Composition and NMDA slope conductance ===\n');
P = base_params(); Bref = mg_block(-55, P);
for fN = [0.25 0.5 0.75]
    r = fN/(1-fN)/Bref; fprintf('fN=%.2f: gN/gA %.2f, peak unblocked conductance ratio %.2f\n', fN, r, r*P.tau_A/P.tau_N);
end
fprintf('B(-55 mV) = %.3f\n', Bref);
for m = [0.5 1 2]
    P.Mg = m; V = linspace(-90, 0, 90001); [B, dB] = mg_block(V, P); sl = B + dB.*(V - P.E_E);
    fprintf('Mg %.1f mM: NSC for V < %.1f mV; slope/chord at -60 mV %.2f\n', m, V(find(sl > 0, 1)), sl(abs(V+60)<1e-6)/B(abs(V+60)<1e-6));
end
% weights after capping
P = base_params(); rs = []; w = [];
for sd = 1:3, C = build_connectivity(P, sd); rs = [rs; sum(C.JE,2); sum(C.JI,2)]; x = C.JE(C.JE>0)*P.KE; w = [w; x]; end %#ok<AGROW>
fprintf('weights after cap: mean %.2f/K, CV %.2f; mean row sum %.2f\n', mean(w), std(w)/mean(w), mean(rs));

fprintf('\n=== Continuation of the quiet branch (Experiment 3b) ===\n');
d = dir(fullfile(root, 'results', 'exp3b_continuation', '*.mat'));
for k = 1:numel(d)
    L = load(fullfile(root, 'results', 'exp3b_continuation', d(k).name));
    j = find(L.maxRe(:,3) > 0, 1); rk = NaN; if ~isempty(j), rk = L.rhoK(j); end
    rkall = arrayfun(@(t) L.rhoK(find(L.maxRe(:,t) > 0, 1)), find(isfinite(L.gE_unstable)));
    fprintf('%s: first unstable gE (tauN 25/50/100/200) %s, restabilised %s, complex %s, rhoK at loss (100 ms) %.3f, rhoK at all Hopf points %s, branch end %.5f (largest real eig of K %.4f, rel. min sing. value %.1e)\n', ...
        d(k).name, mat2str(L.gE_unstable, 4), mat2str(L.gE_restab, 4), mat2str(L.cx_at_loss), rk, mat2str(rkall, 3), L.gE_end, L.realK_end, L.sigmin_end);
end
fprintf('\n--- direct simulation (Exp. 1) versus linear edge (continuation, tau_N = 100 ms, 1 mM) ---\n');
for sd = 1:3, for fN = [0 0.25 0.5 0.75]
    fn = fullfile(root, 'results', 'exp3b_continuation', sprintf('seed%d_fN%02d_Mg100.mat', sd, round(100*fN)));
    if ~exist(fn, 'file'), continue; end
    L = load(fn);
    if isfinite(L.gE_unstable(3)), e = L.gE_unstable(3); elseif isfinite(L.gE_end), e = L.gE_end; else, e = L.gEs(end); end
    s1 = E1.pt(:,3)==sd & E1.pt(:,2)==fN & E1.rateE<4;
    % last coupling of the contiguous quiet range starting at gE = 1.2
    gq = sort(E1.pt(s1,1))'; lastq = NaN; g = 1.2;
    while any(abs(gq - g) < 1e-9), lastq = g; g = g + 0.2; end
    fprintf('seed %d fN %.2f: last quiet gE (direct) %.1f, linear edge %.3f, gap %.0f%%\n', sd, fN, lastq, e, 100*(e - lastq)/e);
end, end
fprintf('\n--- NSC contribution at active-run fixed points of Exp. 3 ---\n');
for m = [1 2]
    s1 = abs(E3.pt(:,2)-m)<1e-9 & E3.has_fp==1 & E3.rateE>=4;
    fprintf('Mg %.0f mM: %s\n', m, mat2str(sort(E3.rhoK(s1)-E3.rhoK_noNSC(s1))', 2));
end
f2 = fullfile(root, 'results', 'validation', 'tangent_map_validation.mat');
if exist(f2, 'file'), L = load(f2); fprintf('\ntangent-map validation [tangent, two-trajectory]: %s\n', mat2str(L.V, 3)); end
end

function sd = f2(E1, E2, E3, k2, k3)
sd = [E1.pt(:,3); E2.pt(k2,3); E3.pt(k3,3)];
end

function [est, ci] = cboot(fun, seed, nb)
% statistic fun(indices) with 95% CI from resampling connectivities
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

function g = g50_of(E1, a1, i, fN)
s = i(E1.pt(i,2) == fN);
X = [ones(numel(s),1) E1.pt(s,1)]; y = double(a1(s)); b = [0;0]; L = diag([0 1e-2]);
for it = 1:100
    p = 1./(1+exp(-X*b)); W = p.*(1-p) + 1e-9;
    db = (X'*(X.*W) + L) \ (X'*(y-p) - L*b); b = b + db; if max(abs(db)) < 1e-10, break; end
end
g = -b(1)/b(2);
end

function r = corr_rank(a, b)
a = a(:); b = b(:); ok = isfinite(a) & isfinite(b); a = a(ok); b = b(ok);
[~, ia] = sort(a); ra(ia) = 1:numel(a); [~, ib] = sort(b); rb(ib) = 1:numel(b);
c = corrcoef(ra, rb); r = c(1,2);
end
