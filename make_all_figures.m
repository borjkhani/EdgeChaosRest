function make_all_figures()
%MAKE_ALL_FIGURES  Regenerate Figures 2-7 and S1-S4 (Figure 1 is drawn in LaTeX/TikZ)
%   from the saved results, and print every number quoted in the text.
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, 'src'), fullfile(here, 'figures_code'));
old = cd(here); c = onCleanup(@() cd(old));
fig2_nmda_biophysics();
fig3_regimes();
fig4_theory();
fig5_nmda_controls();
fig6_hysteresis();
fig7_vision();
figS1_memory_controls();
figS2_homeostasis();
figS3_restoration();
figS4_arclength();
stats_for_text();
stats_controls();
fprintf('Figures written to %s\n', fullfile(here, 'figures'));
end
