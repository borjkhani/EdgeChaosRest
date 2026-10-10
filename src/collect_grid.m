function T = collect_grid(name)
%COLLECT_GRID  Gather results/<name>/pt_*.mat into a struct array T.
d = dir(fullfile('results', name, 'pt_*.mat'));
T = [];
for k = 1:numel(d)
    s = load(fullfile('results', name, d(k).name));
    R = s.R;
    if ~isfield(R, 'ok') || ~R.ok, continue; end
    R = rmfield(R, 'P');
    if isempty(T), T = R; else
        f = fieldnames(T); g = fieldnames(R);
        for q = 1:numel(g), if ~isfield(T, g{q}), [T.(g{q})] = deal([]); end, end
        for q = 1:numel(f), if ~isfield(R, f{q}), R.(f{q}) = []; end, end
        R = orderfields(R, T(1)); T(end+1) = R; %#ok<AGROW>
    end
end
end
