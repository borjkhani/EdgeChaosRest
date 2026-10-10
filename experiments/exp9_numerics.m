function exp9_numerics()
%EXP9_NUMERICS  Numerical and memory controls on selected runs of Exp. 1
%   (connectivities 1-3).
%   Points: the 6 quiet and the 6 active runs whose autonomous exponent is
%   closest to zero, plus 2 deep-quiet (most negative exponent) and 2
%   strongly chaotic (largest exponent) runs.
%   For each point:
%     long : autonomous exponent from 1200 bins (24 s) after 200 discarded,
%            in 4 consecutive blocks of 300 bins (mean and s.d. over blocks)
%     dt   : all measures with the integration step halved (0.25 ms)
%     mem  : memory with 4x longer training and test sets (4800/2400 bins)
%            and the chance baseline for an input that was never presented
%   No-recurrence baseline: the same networks with J^E = J^I = 0, with
%   (a) unchanged drive and (b) background drive raised so that the mean E
%   rate is 2.1 Hz (healthy rate of the V1 network; median quiet rate of
%   Experiment 1: 1.9 Hz).
%   Saves results/exp9_numerics/numerics.mat
outdir = fullfile('results', 'exp9_numerics');
if ~exist(outdir, 'dir'), mkdir(outdir); end
if exist(fullfile(outdir, 'numerics.mat'), 'file'), return; end
T = collect_grid('exp1_regime');
pt = reshape([T.pt], 3, [])'; la = [T.lam_aut]'; rate = [T.rateE]';
ok = pt(:,3) <= 3; q = rate < 4;
iq = find(ok & q); [~, o] = sort(abs(la(iq))); sel = iq(o(1:6));
ia = find(ok & ~q); [~, o] = sort(abs(la(ia))); sel = [sel; ia(o(1:6))];
[~, o] = sort(la(iq)); sel = [sel; iq(o(1:2))];
[~, o] = sort(la(ia), 'descend'); sel = [sel; ia(o(1:2))];
cfg = scan_cfg();
res = struct([]);
for r = 1:numel(sel)
    i = sel(r); P = set_composition(base_params(), pt(i,2), pt(i,1), -55); s = pt(i,3);
    C = build_connectivity(P, s);
    % long autonomous exponent in blocks
    o0 = simulate_net(P, C, struct('T_bins', 200, 'seed', s));
    lb = zeros(4, 1); x = o0.x;
    for b = 1:4
        ob = simulate_net(P, C, struct('T_bins', 300, 'lyap', true, 'lyap_skip', 0, 'seed', s + 10*b, 'x0', x));
        lb(b) = ob.lam; x = ob.x;
    end
    % halved step
    Ph = P; Ph.dt = 0.25; Mh = evaluate_point(Ph, C, s, cfg);
    % longer memory sets + chance baseline
    cm = cfg; cm.do_aut = false; cm.n_train = 4800; cm.n_test = 2400; cm.extra = true;
    Mm = evaluate_point(P, C, s, cm);
    res(r).pt = pt(i,:); res(r).ref = struct('lam_aut', la(i), 'lam_cond', T(i).lam_cond, 'MC', T(i).MC, 'rateE', rate(i));
    res(r).lam_blocks = lb; res(r).dt = Mh; res(r).mem = Mm;
    fprintf('%s: lamA %.2f | blocks %.2f+-%.2f | dt/2 lamA %.2f lamC %.2f (ref %.2f) MC %.2f (ref %.2f) rate %.2f (ref %.2f) | long MC %.2f, chance %.3f\n', ...
        mat2str(pt(i,:), 3), la(i), mean(lb), std(lb), Mh.lam_aut, Mh.lam_cond, T(i).lam_cond, Mh.MC, T(i).MC, Mh.rateE, rate(i), Mm.MC, Mm.MC_shuf);
end
% no-recurrence baseline
base = struct([]);
for s = 1:3
    P = set_composition(base_params(), 0.5, 1.4, -55); C = build_connectivity(P, s);
    C0 = C; C0.JE = 0 * C.JE; C0.JI = 0 * C.JI;
    ce = cfg; ce.do_aut = false; ce.extra = true;
    Ma = evaluate_point(P, C0, s, ce);
    randn('state', 77 + s); lo = -10; hi = 30;   % bisection on I_bg for a mean E rate of 2.1 Hz
    for it = 1:14
        Pb = P; Pb.I_bg = (lo + hi) / 2;
        ob = simulate_net(Pb, C0, struct('T_bins', 300, 'noise', true, 'u', randn(300, 1), 'seed', s));
        if ob.rateE > 2.1, hi = Pb.I_bg; else, lo = Pb.I_bg; end
    end
    Pb = P; Pb.I_bg = (lo + hi) / 2; Mb = evaluate_point(Pb, C0, s, ce);
    base(s).same_drive = Ma; base(s).matched = Mb; base(s).I_bg = Pb.I_bg;
    fprintf('no recurrence, seed %d: same drive MC %.2f NLC %.2f rate %.2f | matched rate (I_bg %.2f) MC %.2f NLC %.2f rate %.2f\n', ...
        s, Ma.MC, Ma.NLC, Ma.rateE, Pb.I_bg, Mb.MC, Mb.NLC, Mb.rateE);
end
save(fullfile(outdir, 'numerics.mat'), 'res', 'base', 'sel');
end
