# NMDA kinetics, Mg²⁺ block and the edge of instability in a biophysical cortical network

[![MATLAB](https://img.shields.io/badge/MATLAB-R2020a%2B-orange.svg)](https://www.mathworks.com/products/matlab.html)
[![GNU Octave](https://img.shields.io/badge/GNU%20Octave-8.x-blue.svg)](https://octave.org)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)

Code and data for the manuscript

> **Separable roles of NMDA receptor kinetics and magnesium block at the edge of instability in a biophysical cortical network, with implications for vision restoration**
> Mehdi Borjkhani, Hadi Borjkhani, Morteza A. Sharif

Every figure and every number in the manuscript is produced by the code in this repository. The saved simulation results are included, so all figures and statistics can be regenerated in minutes without rerunning the simulations.

---

## The model

A conductance-based rate model of a local cortical circuit:

- 160 excitatory and 40 inhibitory units with explicit leak, AMPA, NMDA and GABA_A currents. There are no other currents.
- NMDA currents have the voltage-dependent Mg²⁺ block of Jahr & Stevens (1990).
- Connectivity is sparse (20 excitatory and 5 inhibitory inputs per unit), heterogeneous (lognormal weights) and obeys Dale's law.
- The AMPA/NMDA split is set by an NMDA fraction `f_N` at a fixed total excitatory conductance `gE_tot` at −55 mV.
- Voltages are integrated with exponential Euler (Δt = 0.5 ms).
- Lyapunov exponents use the exact tangent map of that update, in the full 560-dimensional state space.
- Memory is the linear memory capacity of Jaeger (2002) and the quadratic information-processing capacity of Dambre et al. (2012).
- A slow multiplicative factor θ on all excitatory synapses implements homeostatic synaptic scaling. The feedforward drive is scaled by `f_in` to model retinal degeneration (`f_in = 0.25`) and restoration.

The equations are in the header of `src/default_params.m` and in Section 2 of the manuscript.

## Main results

1. **Exact result.** Static (zero-eigenvalue) instabilities of any fixed point are set by a gain matrix *K* that contains no membrane or synaptic time constant. NMDA receptors enter *K* through Mg²⁺-gated coupling and a negative slope conductance.
2. **How the quiet state is lost.**
   - Newton and pseudo-arclength continuation show that without NMDA receptors the quiet state ends at a fold.
   - NMDA receptors can add an earlier oscillatory (Hopf) instability, which NMDA kinetics move by at most 15 %.
   - Raising Mg²⁺ from 0.5 to 2 mM moves the edge about two-fold.
   - NMDA receptors also create a coexisting active state. Bistability of the noise-free network is confirmed up to 20 % below the linear edge.
   - Across ten connectivities, the coupling at which the active state becomes more likely than the quiet one falls by 11–24 % with increasing NMDA fraction (logistic fit with cluster-bootstrap confidence intervals).
3. **Lyapunov exponents and memory dissociate.**
   - Away from the transition, the exponent of the quiet state stays near −c/τ_N (0 < c < 1). It changes four-fold with the NMDA decay time while memory capacity changes by about a quarter.
   - As τ_N grows, the slowest mode's component power decreases relative to the full input-to-rate impulse response (linearised analysis).
4. **Memory is highest in ordered states, not at λ = 0.**
   - Across 520 simulations, linear memory is 3.4 in quiet and 1.2 in active states (difference 2.2, 95 % CI 1.7–2.6), and falls as the conditional exponent approaches zero. Quadratic memory behaves the same way.
   - Memory increases with the distance below each network's own stability edge.
   - Networks without recurrence retain almost no memory.
5. **Visual cortex (exploratory).**
   - After loss of most feedforward drive, homeostatic scaling has no steady state at its set point. With every scaling rule tested, it cycles through hyperactive episodes.
   - Across sampled states of this cycle, restoring input gives hyperactive, chaotic dynamics when θ exceeds a network-specific threshold.
   - Raised Mg²⁺ or partial NMDA block raised that threshold in two of three NMDA-containing networks.

## Quick start

```bash
git clone https://github.com/borjkhani/EdgeChaosRest.git
cd EdgeChaosRest
```

In MATLAB (or GNU Octave), from the repository folder:

```matlab
% regenerate Figs 2-7 and S1-S4 from the saved results, and print every
% number quoted in the manuscript (about 10 min)
make_all_figures

% rerun all simulations (resumable: finished points are skipped)
run_all_experiments
% or share the grid scans between two sessions:
%   session 1: run_all_experiments(1, 2)
%   session 2: run_all_experiments(2, 2)
```

From the command line with Octave:

```bash
octave --no-gui --eval "make_all_figures"
```

Figures are written to `figures/` (PDF and PNG). The numbers quoted in the text are printed by `figures_code/stats_for_text.m` and `figures_code/stats_controls.m`. Their output is saved in `results/stats_for_text_output.txt` and `results/stats_controls_output.txt`.

## Figures and the scripts that make them

| Figure | Content | Script | Experiments |
|---|---|---|---|
| 1 | Model schematic | drawn in LaTeX (TikZ) | – |
| 2 | NMDA biophysics: Mg²⁺ block, negative slope conductance | `fig2_nmda_biophysics.m` | – |
| 3 | Gain × NMDA fraction: phase diagram, quiet vs active states, memory vs λ_cond | `fig3_regimes.m` | 1, 2, 3 |
| 4 | How the quiet state loses stability (continuation, ρ_K, Hopf vs fold) | `fig4_theory.m` | 2, 3, 3b |
| 5 | Mg²⁺ and NMDA kinetics; slowest mode and the input signal | `fig5_nmda_controls.m` | 2, 3, 13 |
| 6 | Steady states of the V1 network versus homeostatic scaling | `fig6_hysteresis.m` | 5 |
| 7 | Degeneration, homeostatic cycling and restoration | `fig7_vision.m` | 6, 11 |
| S1 | Memory controls: distance to edge, delay-resolved and quadratic memory, baselines | `figS1_memory_controls.m` | 1–3, 3b, 9, 10 |
| S2 | Robustness of homeostatic cycling to the scaling rule | `figS2_homeostasis.m` | 8 |
| S3 | Restoration from actual states of the homeostatic cycle | `figS3_restoration.m` | 11 |
| S4 | Pseudo-arclength continuation through the fold | `figS4_arclength.m` | 3c |

## Experiments

| No. | Script | What it does | Approx. time (1 core, Octave) |
|---|---|---|---|
| 1 | `exp1_regime.m` | gain × NMDA fraction, 10 connectivities (400 runs) | 5 h |
| 2 | `exp2_kinetics.m` | gain × NMDA decay time constant (25–200 ms) | 1 h |
| 3 | `exp3_magnesium.m` | gain × [Mg²⁺] (0.25–2 mM) | 1 h |
| 3b | `exp3b_continuation.m` | Newton continuation of the quiet fixed point (fold vs Hopf) | 45 min |
| 3c | `exp3c_arclength.m` | pseudo-arclength continuation through the fold | 15 min |
| 4 | `exp4_tau_heterogeneity.m` | heterogeneous NMDA kinetics (GluN2A/2B mixture) | 20 min |
| 5 | `exp5_hysteresis.m` | V1 steady states versus homeostatic scaling | 40 min |
| 6 | `exp6_vision.m` | degeneration → homeostatic scaling → restoration | 1.5 h |
| 7 | `exp7_current_based.m` | robustness check with current-based AMPA/GABA_A | 20 min |
| 8 | `exp8_homeostasis_robustness.m` | homeostatic rule variants (gain, bounded step, quasi-continuous, floor) | 5 h |
| 9 | `exp9_numerics.m`, `exp9b_replicates.m` | step size, run length, longer readouts, no-recurrence baseline, realisation replicates | 2 h |
| 10 | `exp10_nonlinear.m` | quadratic memory and chance level for all Exp. 1 runs | 2 h |
| 11 | `exp11_restore_states.m` | restoration from stored states of the homeostatic cycle | 1 h |
| 12 | `exp12_coexistence.m` | second attractor below the linear edge (multiple initial states) | 5 min |
| 13 | `exp13_slow_mode.m` | modal analysis of the linearised quiet state | 2 min |
| – | `validate_tangent_map.m` | tangent-map versus two-trajectory Lyapunov exponents | 5 min |

Grid scans (Exp. 1, 2, 3, 4, 7, 10) store one file per point in `results/<experiment>/pt_XXXX.mat` and can be interrupted and resumed.

## Repository structure

```
EdgeChaosRest/
├── run_all_experiments.m     % runs every simulation (resumable; can be split across sessions)
├── make_all_figures.m        % regenerates Figs 2-7, S1-S4 and prints every number quoted in the text
├── src/                      % model and analysis
│   ├── default_params.m      %   all parameters + model equations (header)
│   ├── base_params.m         %   parameter set used in the paper
│   ├── build_connectivity.m  %   sparse lognormal E/I connectivity (Dale's law)
│   ├── simulate_net.m        %   exponential-Euler integration + exact tangent map (Lyapunov)
│   ├── evaluate_point.m      %   autonomous/conditional exponents, memory capacity, baselines
│   ├── memory_capacity.m     %   linear memory capacity (Jaeger 2002), ridge readout
│   ├── nonlinear_capacity.m  %   quadratic information-processing capacity (Dambre et al. 2012)
│   ├── jacobian_fp.m         %   analytic Jacobian at a fixed point
│   ├── newton_fp.m           %   fixed points by Newton's method (stable or unstable)
│   ├── fp_residual.m         %   fixed-point equations and their Jacobian (arclength continuation)
│   ├── gain_criterion.m      %   static gain matrix K, rho_K, NSC and Mg2+ ablations
│   ├── homeostatic_scaling.m %   homeostatic scaling (population-level or cell-autonomous)
│   ├── set_composition.m     %   AMPA/NMDA split at fixed excitatory conductance
│   ├── mg_block.m, rate_fn.m, syn_mode_coeffs.m, relax_to_fp.m
│   └── run_point.m, run_grid.m, collect_grid.m   % scan infrastructure
├── experiments/              % one script per experiment (numbers match the manuscript)
├── figures_code/             % one function per figure, plotting helpers,
│                             %   stats_for_text.m, stats_controls.m, logistic_edge.m
├── figures/                  % Fig2-Fig7, FigS1-FigS4 (.pdf, .png)
├── results/                  % saved results of all experiments (.mat) and statistics output
└── manuscript/               % LaTeX source (Overleaf-ready), references.bib, figures/
```

## Reproducibility notes

- The saved results were produced with GNU Octave 8.4. The same code runs unchanged in MATLAB R2020a or later. Figures exported from MATLAB look sharper than the Octave/gnuplot versions.
- All random numbers (connectivity, input, noise) are seeded, so non-chaotic runs reproduce exactly. In chaotic runs, floating-point differences between machines or sessions (e.g. multithreaded BLAS) can change individual values slightly. When the 400 driven runs of Experiment 1 were repeated in a separate session, 333 were identical; the other 67 (almost all chaotic) changed memory capacity by a median of 0.08, and no run changed state.
- Uncertainties in the manuscript are obtained by resampling connectivities (cluster bootstrap), because runs that share a connectivity are not independent.

## Requirements

MATLAB R2020a or later (no toolboxes needed), or GNU Octave 8.x.

## Citation

If you use this code, please cite:

```bibtex
@unpublished{borjkhani2026nmda,
  title  = {Separable roles of NMDA receptor kinetics and magnesium block at the edge
            of instability in a biophysical cortical network, with implications for
            vision restoration},
  author = {Borjkhani, Mehdi and Borjkhani, Hadi and Sharif, Morteza A.},
  year   = {2026},
  note   = {Manuscript submitted},
  url    = {https://github.com/borjkhani/EdgeChaosRest}
}
```

## Contact

Mehdi Borjkhani, International Centre for Translational Eye Research (ICTER), Institute of Physical Chemistry, Polish Academy of Sciences, Warsaw, Poland. mborjkhani@ichf.edu.pl

## License

MIT, see [LICENSE](LICENSE).
