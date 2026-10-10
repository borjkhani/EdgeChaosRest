function S = fig_style()
%FIG_STYLE  Shared colours, markers and sizes for all manuscript figures.
%   Categorical colours are used in a fixed order (slot 1..4) and always
%   paired with a distinct marker so that identity never relies on colour.
S.col = [0.165 0.471 0.839;   % 1 blue   #2a78d6
         0.922 0.408 0.204;   % 2 orange #eb6834
         0.106 0.686 0.478;   % 3 aqua   #1baf7a
         0.929 0.631 0.000;   % 4 yellow #eda100
         0.384 0.314 0.839];  % 5 violet #6250d6 (reference lines only)
S.mk  = {'o', 's', '^', 'd', 'v'};
S.ink = [0.04 0.04 0.04]; S.ink2 = [0.32 0.32 0.31]; S.grid = [0.85 0.85 0.85];
S.lw = 1.5; S.ms = 5; S.fs = 9;
end
