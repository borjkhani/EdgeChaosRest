function D = load_scan(name)
%LOAD_SCAN  Results of a parameter scan as column arrays.
%   D.pt (n x 3: swept value 1, swept value 2, seed) and one column per
%   scalar measure (lam_aut, lam_cond, MC, rhoK, ...).
here = fileparts(mfilename('fullpath'));
old = cd(fullfile(here, '..'));
addpath(fullfile(here, '..', 'src'));
T = collect_grid(name);
cd(old);
D.pt = reshape([T.pt], numel(T(1).pt), [])';
f = {'lam_aut','lam_cond','MC','rhoK','rhoK_noNSC','rhoK_noMg','maxReJ','rateE', ...
     'sd_aut','sync','t_mem50','has_fp','NSC','G','V'};
for q = 1:numel(f)
    v = nan(numel(T), 1);
    for k = 1:numel(T)
        x = T(k).(f{q}); if ~isempty(x), v(k) = double(x); end
    end
    D.(f{q}) = v;
end
D.MCk = reshape([T.MCk], [], numel(T))';
end
