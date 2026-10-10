function run_grid(name, pts, make_P, cfg, worker, nworkers)
%RUN_GRID  Evaluate a list of grid points, sharing them between workers.
%   name     : results/<name>/ receives one pt_XXXX.mat per point
%   pts      : n x d matrix; column d is the connectivity seed
%   make_P   : function handle, P = make_P(pts(r,:))
%   cfg      : options passed to run_point / evaluate_point
%   worker, nworkers : this process handles rows worker, worker+nworkers, ...
%   Points whose result file already exists are skipped, so a scan can be
%   resumed or split across MATLAB sessions.
if nargin < 5, worker = 1; nworkers = 1; end
outdir = fullfile('results', name);
if ~exist(outdir, 'dir'), mkdir(outdir); end
for r = worker:nworkers:size(pts, 1)
    fn = fullfile(outdir, sprintf('pt_%04d.mat', r));
    if exist(fn, 'file'), continue; end
    P = make_P(pts(r, :));
    t0 = tic;
    try
        R = run_point(P, pts(r, end), cfg);
        R.ok = true;
    catch err
        R = struct('ok', false, 'msg', err.message);
    end
    R.pt = pts(r, :); R.P = P; R.time = toc(t0);
    save(fn, 'R');
    if R.ok
        fprintf('[%s %d/%d] %s  lamA=%6.2f lamC=%6.2f MC=%5.2f rhoK=%5.2f (%.0fs)\n', name, r, size(pts,1), ...
            mat2str(pts(r,:), 3), R.lam_aut, R.lam_cond, R.MC, R.rhoK, R.time);
    else
        fprintf('[%s %d] FAILED: %s\n', name, r, R.msg);
    end
end
end
