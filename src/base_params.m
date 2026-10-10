function P = base_params()
%BASE_PARAMS  Parameter set used for all figures of the manuscript.
%   Starts from default_params and sets the sparse, heterogeneous E/I
%   connectivity and input statistics used in the paper. Every experiment
%   modifies only the parameters it sweeps.
P = default_params();
P.kappa   = 3;      % mV, rate-function slope (expansive regime at low rates)
P.KE      = 20;     % excitatory in-degree
P.KI      = 5;      % inhibitory in-degree
P.cv_w    = 2;      % lognormal weight CV (capped at 5x mean)
P.tau_G   = 5;      % ms
P.gG      = 20;     % GABA_A conductance (x g_L)
P.I_bg    = 8;      % mV, background drive
P.I_ff    = 6;      % mV, mean feedforward (visual) drive at f_in = 1
P.sig_u   = 1.5;    % mV, input signal amplitude
P.sig_eta = 0.3;    % mV, private noise
P.V_ref   = -55;    % mV, reference voltage for AMPA/NMDA composition
end
