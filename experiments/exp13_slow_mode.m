function exp13_slow_mode()
%EXP13_SLOW_MODE  Why does the slowest mode set the Lyapunov exponent but not
%   memory? At every quiet fixed point of Experiment 2 (fN = 0.5,
%   tau_N = 25-200 ms, gE_tot = 1.2 and 1.4, connectivities 1-2) the
%   network is linearised, dx/dt = J x + b xi(t), with b the input direction
%   (sig_u w / tau_m in the voltage components, 0 elsewhere) and xi white
%   noise (a linear, white-input approximation of the piecewise-constant
%   input of the memory task). The read-out are the rates, y = C x with
%   C = [Phi' 0 0 0].
%   With J = V diag(mu) V^-1 (any eigenvector scaling; the results do not
%   depend on it), the stationary covariance is Sigma = V S V^H with
%       S_ij = -(l_i b) conj(l_j b) / (mu_i + conj(mu_j)),  l_i = rows of V^-1,
%   which includes all cross terms between non-orthogonal modes, and the
%   lagged rate covariance is C V diag(exp(mu tau)) S V^H C'.
%   The memory-relevant quantity is the input signal carried by the rates at
%   lag tau, i.e. the impulse response h(tau) = C exp(J tau) b
%   = C V diag(exp(mu tau)) V^-1 b. Reported for the slowest mode k1 (and its
%   complex conjugate, if any):
%     decay   : Re mu_k1 (1/s)
%     imp     : share of the impulse-response power sum_tau ||h(tau)||^2,
%               tau = 20, 40, ..., 1000 ms, carried by the slow mode,
%               ||h_k1(tau)||^2 summed / ||h(tau)||^2 summed (cross terms
%               between modes are part of ||h||^2)
%     imp_lag : the same at single lags of 20, 60, 100 and 200 ms
%     share0  : change of the total stationary rate variance trace(C Sigma C')
%               when the input to the slow mode is removed, relative to the
%               total; because modes are not orthogonal this is not an
%               additive share and can be negative
%     iso     : isolated-mode measure |l_k1 b|^2 ||C v_k1||^2 / (2|Re mu_k1|)
%               relative to its sum over modes (ignores cross terms; for
%               comparison only)
%     resid   : relative residual of the Lyapunov equation (check)
%   Saves results/exp13_slow_mode/slow_mode.mat
outdir = fullfile('results', 'exp13_slow_mode');
if ~exist(outdir, 'dir'), mkdir(outdir); end
gEs = [1.2 1.4]; tauN = [25 50 100 200]; out = []; curves = {};
for s = 1:2, for g = gEs, for t = tauN
    P = set_composition(base_params(), 0.5, g, -55); P.tau_N = t;
    C = build_connectivity(P, s);
    [x, conv] = relax_to_fp(P, C, [], 250);
    if conv > 1e-4 || P.r_max * mean(x.a) > 4, continue; end
    J = jacobian_fp(P, C, x) * 1000;                 % 1/s
    [Vv, D] = eig(J); mu = diag(D); L = inv(Vv); n = size(J, 1); N = C.N;
    b = zeros(n, 1); b(1:N) = 1000 * P.sig_u * C.w / P.tau_m;
    [~, dnu] = rate_fn(x.V, P); Cm = [diag(dnu), zeros(N, n - N)];
    bt = L * b; CV = Cm * Vv;
    [~, k1] = max(real(mu)); ks = find(abs(mu - conj(mu(k1))) < 1e-9 | (1:n)' == k1);
    tot = @(btv, tau) real(trace(CV * diag(exp(mu * tau)) * (-(btv * btv') ./ (mu + mu')) * CV'));
    bw = bt; bw(ks) = 0;
    sh0 = 1 - tot(bw, 0) / tot(bt, 0);
    h = @(btv, tau) CV * (exp(mu * tau) .* btv);
    T = 0.02:0.02:1; num = 0; den = 0;
    for tau = T, hv = h(bt, tau); hs = hv - h(bw, tau); num = num + norm(hs)^2; den = den + norm(hv)^2; end
    lagr = arrayfun(@(tau) norm(h(bt, tau) - h(bw, tau))^2 / norm(h(bt, tau))^2, [0.02 0.06 0.1 0.2]);
    curve = arrayfun(@(tau) norm(h(bt, tau))^2, T); curve_wo = arrayfun(@(tau) norm(h(bw, tau))^2, T);
    iso = abs(bt).^2 .* sum(abs(CV).^2, 1)' ./ (2 * abs(real(mu)));
    Sig = Vv * (-(bt * bt') ./ (mu + mu')) * Vv';
    res = norm(J * Sig + Sig * J' + b * b', 'fro') / norm(b * b', 'fro');
    out(end+1, :) = [s g t real(mu(k1)) num/den lagr sh0 sum(iso(ks)) / sum(iso) res]; %#ok<AGROW>
    curves{size(out, 1)} = [curve; curve_wo]; %#ok<AGROW>
    fprintf('seed %d gE %.1f tauN %3d: slowest mode %.2f /s; slow-mode share of impulse-response power %.3f (lags 20/60/100/200 ms: %.3f %.3f %.3f %.3f); change of stationary variance %.3f; isolated-mode measure %.1e; Lyapunov residual %.1e\n', out(end, :));
end, end, end
lags = 0.02:0.02:1;
save(fullfile(outdir, 'slow_mode.mat'), 'out', 'curves', 'lags');
end
