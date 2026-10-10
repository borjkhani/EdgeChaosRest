function cfg = scan_cfg()
%SCAN_CFG  Simulation lengths shared by all parameter scans.
cfg = struct('T_aut', 300, 'skip_aut', 100, 'n_wash', 200, ...
             'n_train', 1200, 'n_test', 600, 'K', 50);
end
