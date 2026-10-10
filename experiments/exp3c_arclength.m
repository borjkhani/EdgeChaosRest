function exp3c_arclength(seed, fN)
%EXP3C_ARCLENGTH  Pseudo-arclength continuation of the quiet fixed point in
%   gE_tot (1 mM Mg2+), to check that the end of the Newton branch
%   (Experiment 3b) is a fold: the branch must reach a maximum of gE_tot
%   and turn back as an unstable (saddle) branch. Unknowns y = (V, gE_tot);
%   predictor along the tangent of F(V, g) = 0, Newton corrector on
%   [F; t'(y - y_pred)] = 0, with g scaled by 100 in the arclength norm.
%   Along the branch: gE_tot, mean V, mean E rate, smallest singular value
%   of dF/dV, largest real eigenvalue of K (real eigenvalues only) and the
%   number of Jacobian eigenvalues with positive real part (tau_N = 100 ms).
%   Saves results/exp3c_arclength/seed<s>_fN<100 fN>.mat
outdir = fullfile('results', 'exp3c_arclength');
if ~exist(outdir, 'dir'), mkdir(outdir); end
P0 = base_params(); C = build_connectivity(P0, seed);
Pf = @(g) set_composition(P0, fN, g, -55);
sg = 100; N = C.N;
g = 1.0; [x, ok] = newton_fp(Pf(g), C, (P0.E_L + 8) * ones(N, 1)); V = x.V;
assert(ok);
tv = tangent(Pf, C, V, g, sg, []);
ds = 0.5; out = [];
for step = 1:4000
    yp = [V; sg * g] + ds * tv;
    y = yp; conv = false;
    for it = 1:30
        Vc = y(1:N); gc = y(end) / sg;
        [F, JF] = fp_residual(Pf(gc), C, Vc);
        Fg = (fp_residual(Pf(gc + 1e-6), C, Vc) - fp_residual(Pf(gc - 1e-6), C, Vc)) / 2e-6;
        Jx = [JF, Fg / sg; tv'];
        r = [F; tv' * (y - yp)];
        dy = -Jx \ r; y = y + dy;
        if max(abs(dy)) < 1e-9 && max(abs(F)) < 1e-8, conv = true; break; end
    end
    if ~conv || any(~isfinite(y))
        ds = ds / 2; if ds < 1e-4, break; end, continue
    end
    V = y(1:N); g = y(end) / sg;
    tv = tangent(Pf, C, V, g, sg, tv);
    P = Pf(g); [~, JF] = fp_residual(P, C, V);
    xs = struct('V', V, 'a', rate_fn(V(C.isE), P), 'n', rate_fn(V(C.isE), P), 'b', rate_fn(V(~C.isE), P));
    [Jac, info] = jacobian_fp(P, C, xs); ev = eig(Jac);
    D = diag(info.dnu); B = mg_block(V, P);
    rowE = P.gA * info.dfE + P.gN * B .* info.dfN;
    K = [diag(rowE ./ info.Gloc) * C.JE * D(C.isE, C.isE), diag(P.gG * info.dfI ./ info.Gloc) * C.JI * D(~C.isE, ~C.isE)];
    ek = eig(K); er = real(ek(abs(imag(ek)) < 1e-9));
    out(end+1, :) = [g, mean(V), P.r_max * mean(xs.a), min(svd(JF)) / max(svd(JF)), max(er), sum(real(ev) > 0)]; %#ok<AGROW>
    if size(out, 1) > 5 && (g < max(out(:,1)) - 0.3 || out(end, 3) > 30 || g > 3.6), break; end
    ds = min(2, ds * 1.2);
end
% first turning point of gE_tot (fold) and stability around it
dg = diff(out(:,1)); i1 = find(dg < 0, 1);   % first point after which gE decreases
g_fold = out(i1, 1);
i_unst = find(out(:,6) > 0, 1); g_loss = out(max(i_unst - 1, 1), 1);
% stable low-rate segments anywhere on the computed branch
st = out(:,6) == 0 & out(:,3) < 4; g_lastq = max(out(st, 1));
gmax = max(out(:,1));
fprintf('seed %d fN %.2f: first loss of stability at gE %.4f; first fold at gE %.4f (rel. sing. %.1e, real eig K %.3f), unstable eigenvalues after fold %d; largest gE of any stable low-rate (<4 Hz) fixed point %.4f\n', ...
    seed, fN, g_loss, g_fold, out(i1,4), out(i1,5), out(min(i1+1, end),6), g_lastq);
save(fullfile(outdir, sprintf('seed%d_fN%02d.mat', seed, round(100*fN))), 'out', 'gmax', 'g_fold', 'g_loss', 'g_lastq');
end

function t = tangent(Pf, C, V, g, sg, tprev)
[~, JF] = fp_residual(Pf(g), C, V);
Fg = (fp_residual(Pf(g + 1e-6), C, V) - fp_residual(Pf(g - 1e-6), C, V)) / 2e-6;
z = -JF \ (Fg / sg);            % dV / d(sg g)
t = [z; 1]; t = t / norm(t);
if isempty(tprev)
    return
end
% at a fold z diverges; use the null vector of [JF, Fg/sg] instead
A = [JF, Fg / sg]; [~, ~, W] = svd(A); t = W(:, end);
if t' * tprev < 0, t = -t; end
end
