function exp12_coexistence()
%EXP12_COEXISTENCE  Is the active state seen below the linear edge a second
%   attractor of the noise-free network with constant input?
%   For connectivities 1-3 and every coupling at which the Experiment-1 run
%   was active although the quiet fixed point was still linearly stable
%   (coupling below the edge found by continuation, tau_N = 100 ms, 1 mM),
%   the network is simulated without noise and with constant input (u = 0)
%   from (a) the quiet fixed point (followed by Newton continuation from
%   gE_tot = 1) displaced by a random perturbation of 0.5 mV s.d., (b) the
%   state reached after 2 s of the noisy, input-driven activity of
%   Experiment 1, and (c) the default initial state of the simulations
%   (all voltages at E_L + 15 mV). Final mean E rate and the largest voltage
%   change over the last bin (fixed point if < 1e-4 mV) are recorded.
%   Saves results/exp12_coexistence/coexistence.mat
outdir = fullfile('results', 'exp12_coexistence');
if ~exist(outdir, 'dir'), mkdir(outdir); end
here = pwd; T = collect_grid('exp1_regime');
pt = reshape([T.pt], 3, [])'; rate = [T.rateE]';
rows = [];
for s = 1:3, for f = [0 0.25 0.5 0.75]
    L = load(fullfile(here, 'results', 'exp3b_continuation', sprintf('seed%d_fN%02d_Mg100.mat', s, round(100*f))));
    if isfinite(L.gE_unstable(3)), e = L.gE_unstable(3); else, e = L.gE_end; end
    sel = find(pt(:,3) == s & pt(:,2) == f & rate >= 4 & pt(:,1) < e);
    for i = sel', rows(end+1, :) = [pt(i,:) e]; end %#ok<AGROW>
end, end
out = nan(size(rows, 1), 6);
for r = 1:size(rows, 1)
    P0 = base_params(); C = build_connectivity(P0, rows(r,3));
    V = (P0.E_L + 8) * ones(C.N, 1);
    for g = [1.0:0.005:rows(r,1) rows(r,1)]
        [xq, ok] = newton_fp(set_composition(P0, rows(r,2), g, -55), C, V); V = xq.V;
    end
    P = set_composition(P0, rows(r,2), rows(r,1), -55);
    randn('state', 9 + r); xp = xq; xp.V = xq.V + 0.5 * randn(C.N, 1);
    oq = simulate_net(P, C, struct('T_bins', 400, 'seed', 7, 'x0', xp));
    randn('state', 1000 + rows(r,3));
    od = simulate_net(P, C, struct('u', randn(100, 1), 'noise', true, 'seed', 2000 + rows(r,3)));
    ob = simulate_net(P, C, struct('T_bins', 400, 'seed', 7, 'x0', od.x));
    oa = simulate_net(P, C, struct('T_bins', 400, 'seed', 7));
    out(r, :) = [mean_last(oq, P, C), conv_last(oq, P), mean_last(ob, P, C), conv_last(ob, P), mean_last(oa, P, C), conv_last(oa, P)];
    fprintf('seed %d fN %.2f gE %.1f (edge %.3f; Newton ok %d): near quiet FP %.2f Hz (dV %.1e) | from driven active state %.2f Hz (dV %.1e) | default start %.2f Hz (dV %.1e)\n', ...
        rows(r,3), rows(r,2), rows(r,1), rows(r,4), ok, out(r,:));
end
save(fullfile(outdir, 'coexistence.mat'), 'rows', 'out');
end

function r = mean_last(o, P, C)
r = P.r_max * mean(mean(o.R(C.isE, end-49:end)));
end

function d = conv_last(o, P)
R = o.R(:, end-1:end); V = P.kappa * log(R ./ (1 - R));
d = max(abs(V(:,2) - V(:,1)));
end
