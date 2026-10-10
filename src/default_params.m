function P = default_params()
%DEFAULT_PARAMS  Parameters of the conductance-based E/I rate network.
%
%   tau_m dV_i/dt = -(V_i - E_L)
%                   - gA * (J^E a)_i        * (V_i - E_E)          (AMPA)
%                   - gN * (J^E n)_i * B(V_i) * (V_i - E_E)        (NMDA)
%                   - gG * (J^I b)_i        * (V_i - E_I)          (GABA_A)
%                   + I_i(t)
%   tau_A da_j/dt = -a_j + nu_j,  tau_N dn_j/dt = -n_j + nu_j,
%   tau_G db_j/dt = -b_j + nu_j,  nu_j = phi(V_j)
%
%   Conductances gA, gN, gG are in units of the leak conductance, so the
%   equations are dimensionless except for voltages (mV) and times (ms).
%   I_i(t) is an input expressed as a voltage (input current x leak resistance).
%
%   All parameters are listed here so that every figure can be regenerated
%   from this single file plus the experiment scripts.

% ---------------- network size and connectivity ----------------
P.NE   = 160;          % excitatory units
P.NI   = 40;           % inhibitory units
P.KE   = 40;           % excitatory in-degree (each unit)
P.KI   = 10;           % inhibitory in-degree (each unit)
P.cv_w = 1.0;          % coefficient of variation of (lognormal) weights

% ---------------- single-unit properties ----------------
P.tau_m = 10;          % ms, membrane (population) time constant
P.E_L   = -70;         % mV, leak reversal
P.E_E   = 0;           % mV, AMPA/NMDA reversal
P.E_I   = -80;         % mV, GABA_A reversal
P.theta = -50;         % mV, half-activation of the rate function
P.kappa = 4;           % mV, slope of the rate function
P.r_max = 50;          % Hz, maximal rate (used only to report rates in Hz)

% ---------------- synapses ----------------
P.tau_A = 5;           % ms, AMPA decay
P.tau_N = 100;         % ms, NMDA decay; scalar or NE x 1 vector (per presynaptic E unit)
P.tau_G = 10;          % ms, GABA_A decay
P.gA    = 3.0;         % AMPA conductance (x leak), recurrent E->all
P.gN    = 0.0;         % NMDA conductance (x leak), unblocked; set by experiments
P.gG    = 6.0;         % GABA_A conductance (x leak), recurrent I->all
P.Mg    = 1.0;         % mM, extracellular Mg2+
P.mg_gamma = 3.57;     % mM   (Jahr & Stevens 1990)
P.mg_delta = 0.062;    % 1/mV (Jahr & Stevens 1990)

% ---------------- inputs ----------------
P.I_bg   = 14;         % mV, background (non-visual) drive, all units
P.I_ff   = 6;          % mV, mean feedforward (thalamic/visual) drive
P.sig_u  = 4;          % mV, amplitude of the input signal u(t)
P.sig_eta = 1.0;       % mV, private noise amplitude (white, see simulate_net)
P.f_in   = 1.0;        % feedforward fraction (1 = intact retina/LGN input)

% ---------------- numerics ----------------
P.dt     = 0.5;        % ms, Euler step
P.bin_ms = 20;         % ms, input sample duration and readout bin
end
