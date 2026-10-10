function V = validate_tangent_map()
%VALIDATE_TANGENT_MAP  Compare the tangent-map Lyapunov exponent with an
%   independent two-trajectory (finite-difference) estimate.
%   For three networks (fN = 0.5; gE_tot = 1.2 ordered, 1.6 near the edge,
%   2.4 active; connectivity seed 1) the autonomous exponent is computed
%   (i) with the exact tangent map of simulate_net and (ii) by following a
%   second trajectory displaced by 1e-7 and renormalised every input bin.
%   Both estimates use bins 51-250 after a 100-bin transient (the first 50
%   bins let each perturbation align with the leading direction).
P0 = base_params(); gEs = [1.2 1.6 2.4];
V = zeros(numel(gEs), 2);
for k = 1:numel(gEs)
    P = set_composition(P0, 0.5, gEs(k), -55); C = build_connectivity(P, 1);
    o = simulate_net(P, C, struct('T_bins', 100, 'seed', 2)); x0 = o.x;
    ot = simulate_net(P, C, struct('T_bins', 250, 'lyap', true, 'lyap_skip', 50, 'x0', x0, 'seed', 5));
    d0 = 1e-7; randn('state', 11);
    xa = x0; xb = x0; xb.V = xb.V + d0 * randn(C.N, 1) / sqrt(C.N);
    nb = sqrt(sum((xb.V - xa.V).^2)); s = d0 / nb; xb.V = xa.V + (xb.V - xa.V) * s;
    lsum = 0;
    for t = 1:250
        oa = simulate_net(P, C, struct('T_bins', 1, 'x0', xa)); ob = simulate_net(P, C, struct('T_bins', 1, 'x0', xb));
        xa = oa.x; xb = ob.x;
        d = sqrt(sum((xa.V-xb.V).^2) + sum((xa.a-xb.a).^2) + sum((xa.n-xb.n).^2) + sum((xa.b-xb.b).^2));
        if t > 50, lsum = lsum + log(d / d0); end
        f = d0 / d;
        xb.V = xa.V + (xb.V - xa.V) * f; xb.a = xa.a + (xb.a - xa.a) * f;
        xb.n = xa.n + (xb.n - xa.n) * f; xb.b = xa.b + (xb.b - xa.b) * f;
    end
    V(k, :) = [ot.lam, lsum / (200 * P.bin_ms * 1e-3)];
    fprintf('gE = %.1f: tangent map %.3f 1/s, two trajectories %.3f 1/s\n', gEs(k), V(k, 1), V(k, 2));
end
outdir = fullfile('results', 'validation'); if ~exist(outdir, 'dir'), mkdir(outdir); end
save(fullfile(outdir, 'tangent_map_validation.mat'), 'V', 'gEs');
end
