function run_all_experiments(worker, nworkers)
%RUN_ALL_EXPERIMENTS  Run every simulation of the manuscript.
%
%   run_all_experiments            % everything in one MATLAB session
%   run_all_experiments(w, n)      % share the grid scans between n sessions
%
%   Results go to results/<experiment>/. Points already computed are
%   skipped, so the function can be interrupted and resumed. Approximate
%   single-core times (Octave 8.4; MATLAB is usually faster):
%     exp1_regime            400 points   ~5 h
%     exp2_kinetics           80 points   ~60 min
%     exp3_magnesium          80 points   ~60 min
%     exp3b_continuation      18 branches ~45 min (Newton continuation)
%     exp4_tau_heterogeneity  24 points   ~20 min
%     exp5_hysteresis          8 runs     ~40 min
%     exp6_vision              6 runs     ~90 min
%     exp7_current_based      24 points   ~20 min
%     exp3c_arclength         12 branches ~15 min (pseudo-arclength continuation)
%     exp8_homeostasis       36 runs     ~5 h   (homeostatic rule variants)
%     exp9_numerics          16 points   ~1.5 h (step size, run length, baselines)
%     exp10_nonlinear       400 points   ~2 h   (quadratic memory, chance level)
%     exp11_restore_states    6 runs     ~1 h   (restoration from cycle states)
%     exp12_coexistence      11 points   ~5 min
%     exp13_slow_mode        16 points   ~2 min
%   File numbers match the experiment numbers in the manuscript.
if nargin < 1, worker = 1; nworkers = 1; end
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, 'src'), fullfile(here, 'experiments'));
old = cd(here); c = onCleanup(@() cd(old));
exp1_regime(worker, nworkers);
exp2_kinetics(worker, nworkers);
exp3_magnesium(worker, nworkers);
exp4_tau_heterogeneity(worker, nworkers);
exp7_current_based(worker, nworkers);
if worker == 1
    for s = 1:3
        for f = [0 0.25 0.5 0.75], exp3b_continuation(s, f, 1); end
        exp3b_continuation(s, 0.5, 0.5); exp3b_continuation(s, 0.5, 2);
    end
    validate_tangent_map();
    for cnd = 1:4, for s = 1:2, exp5_hysteresis(cnd, s); end, end
    for s = 1:3, for cnd = [1 2], exp6_vision(s, cnd); end, end
    for s = 1:3, for f = [0 0.25 0.5 0.75], exp3c_arclength(s, f); end, end
    v = {'eta0.3', 'eta0.1', 'eta0.03', 'tanh', 'cont', 'floor'};
    for s = 1:3, for cnd = [1 2], for k = 1:numel(v), exp8_homeostasis_robustness(s, cnd, v{k}); end, end, end
    exp9_numerics(); exp9b_replicates();
    for s = 1:3, for cnd = [1 2], exp11_restore_states(s, cnd); end, end
    exp12_coexistence(); exp13_slow_mode();
end
exp10_nonlinear(worker, nworkers);
end
