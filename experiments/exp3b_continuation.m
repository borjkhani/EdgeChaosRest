function exp3b_continuation(seed, fN, Mg)
%EXP3B_CONTINUATION  Follow the quiet fixed-point branch in gE_tot by Newton
%   continuation and determine how it loses stability, for four NMDA decay
%   constants (Experiment 3b of the manuscript).
%
%   The branch (V*, rho_K) does not depend on tau_N (Proposition 1), so it is
%   computed once; for each tau_N the full Jacobian spectrum at every branch
%   point gives its stability. The branch is followed from gE_tot = 1.0 in
%   steps of 0.005 until Newton fails or the solution jumps by more than
%   3 mV (another branch), or gE_tot = 3.4. The end of the branch is then
%   refined by halving the step down to 1e-5; at the refined end the largest
%   real eigenvalue of K and the smallest singular value of D - M*Phi' are
%   recorded (both indicate a fold: 1 and 0, respectively). Mg: extracellular Mg2+ (mM); conductances are set by
%   set_composition at 1 mM, as in Experiment 3.
%   Saves results/exp3b_continuation/seed<s>_fN<100 fN>_Mg<100 Mg>.mat
if nargin < 1, seed = 1; end
if nargin < 2, fN = 0.5; end
if nargin < 3, Mg = 1; end
outdir = fullfile('results', 'exp3b_continuation');
if ~exist(outdir, 'dir'), mkdir(outdir); end
tauN = [25 50 100 200];
gEs = 1.0:0.005:3.4;
P0 = base_params(); C = build_connectivity(P0, seed);
V = (P0.E_L + 8) * ones(C.N, 1);           % start on the low (quiet) branch
n = numel(gEs);
rhoK = nan(n, 1); rate = nan(n, 1); Vm = nan(n, 1);
maxRe = nan(n, numel(tauN)); isCx = false(n, numel(tauN)); imagPart = nan(n, numel(tauN));
last = 0;
for k = 1:n
    P = set_composition(P0, fN, gEs(k), -55); P.Mg = Mg;   % conductances set at 1 mM
    [x, ok] = newton_fp(P, C, V);
    if ~ok, break; end
    % reject jumps to another branch: the continuation step must be small
    if k > 1 && max(abs(x.V - V)) > 3, break; end
    V = x.V; last = k;
    S = gain_criterion(P, C, x);
    rhoK(k) = S.rhoK; rate(k) = P.r_max * mean(x.a); Vm(k) = mean(x.V);
    for t = 1:numel(tauN)
        Pt = P; Pt.tau_N = tauN(t);
        ev = eig(jacobian_fp(Pt, C, x));
        [~, im] = max(real(ev));
        maxRe(k, t) = 1000 * real(ev(im)); imagPart(k, t) = 1000 * abs(imag(ev(im)));
        isCx(k, t) = abs(imag(ev(im))) > 1e-9;
    end
end
gE_fold = NaN; realK_end = NaN; sigmin_end = NaN; gE_end = NaN;
if last >= 1 && last < n
    % refine the end of the branch by step halving
    g = gEs(last); h = 0.0025; Vr = V;
    while h >= 1e-5
        P = set_composition(P0, fN, g + h, -55); P.Mg = Mg;
        [x2, ok2] = newton_fp(P, C, Vr);
        if ok2 && max(abs(x2.V - Vr)) < 1, g = g + h; Vr = x2.V; else, h = h / 2; end
    end
    P = set_composition(P0, fN, g, -55); P.Mg = Mg;
    [xr, okr] = newton_fp(P, C, Vr);
    if okr
        [~, info] = jacobian_fp(P, C, xr); Bv = info.B; isE = C.isE;
        rowE = P.gA * info.dfE + P.gN * Bv .* info.dfN; rowI = P.gG * info.dfI;
        MP = [diag(rowE) * C.JE * diag(info.dnu(isE)), diag(rowI) * C.JI * diag(info.dnu(~isE))];
        ev = eig(diag(1 ./ info.Gloc) * MP); evr = real(ev(abs(imag(ev)) < 1e-9));
        realK_end = max(evr); sigmin_end = min(svd(diag(info.Gloc) - MP)) / max(info.Gloc);
        gE_end = g; gE_fold = g;
    end
end
gE_unstable = nan(1, numel(tauN)); cx_at_loss = false(1, numel(tauN)); gE_restab = nan(1, numel(tauN));
for t = 1:numel(tauN)
    j = find(maxRe(1:last, t) > 0, 1);
    if ~isempty(j)
        gE_unstable(t) = gEs(j); cx_at_loss(t) = isCx(j, t);
        r = find(maxRe(j:last, t) <= 0, 1);
        if ~isempty(r), gE_restab(t) = gEs(j + r - 1); end
    end
end
fprintf('seed %d fN %.2f Mg %.2f: branch end gE=%.5f (largest real eig of K %.4f, rel. min singular value %.2e); first unstable gE per tauN: %s; complex at loss: %s\n', ...
    seed, fN, Mg, gE_end, realK_end, sigmin_end, mat2str(gE_unstable, 4), mat2str(cx_at_loss));
gEs = gEs(1:last); rhoK = rhoK(1:last); rate = rate(1:last); Vm = Vm(1:last);
maxRe = maxRe(1:last, :); isCx = isCx(1:last, :); imagPart = imagPart(1:last, :);
save(fullfile(outdir, sprintf('seed%d_fN%02d_Mg%03d.mat', seed, round(100 * fN), round(100 * Mg))), ...
    'gEs', 'rhoK', 'rate', 'Vm', 'maxRe', 'isCx', 'imagPart', 'tauN', 'gE_fold', ...
    'gE_unstable', 'cx_at_loss', 'seed', 'fN', 'Mg', 'gE_end', 'realK_end', 'sigmin_end', 'gE_restab');
end
