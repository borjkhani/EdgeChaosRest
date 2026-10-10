function out = simulate_net(P, C, opt)
%SIMULATE_NET  Integrate the conductance-based E/I rate network.
%   Voltages use exponential Euler (exact for conductances frozen over one
%   step, unconditionally stable); gating variables use exact exponential
%   updates. Step size P.dt.
%
%   out = simulate_net(P, C, opt)
%
%   opt fields (all optional):
%     u        : T_bins x 1 input sequence (held constant for P.bin_ms each).
%                Default: zeros(opt.T_bins,1).
%     T_bins   : number of bins (if u not given)
%     noise    : true/false, private white noise of amplitude P.sig_eta
%     lyap     : true/false, compute the maximal (conditional) Lyapunov
%                exponent along the realised trajectory (Benettin method,
%                tangent dynamics of the Euler map, renormalised every bin)
%     lyap_skip: bins discarded before accumulating the exponent
%     x0       : initial state struct with fields V,a,n,b
%     V_init   : mean initial voltage (mV) if x0 is not given (default E_L+15);
%                gating variables start at their steady state for that V
%     seed     : seed for the noise and the initial perturbation
%     scaleE   : N x 1 multiplicative scaling of excitatory (AMPA+NMDA)
%                conductances onto each unit (homeostatic synaptic scaling);
%                default ones.
%
%   out fields:
%     R      : N x T_bins normalised rates sampled at the end of each bin
%     lam    : maximal Lyapunov exponent (1/s); NaN if not requested
%     x      : final state
%     rateE, rateI : mean rates (Hz) of E and I units over the run
%
%   The model is defined in default_params.m. The tangent vector obeys the
%   exact Jacobian of the discrete update used here, so lam is the
%   exponent of the simulated system (it converges to that of the ODE as
%   dt -> 0).

if ~isfield(opt, 'u') || isempty(opt.u)
    opt.u = zeros(opt.T_bins, 1);
end
u      = opt.u(:);
T_bins = numel(u);
noise  = isfield(opt, 'noise') && opt.noise;
dolyap = isfield(opt, 'lyap') && opt.lyap;
skip   = 0; if isfield(opt, 'lyap_skip'), skip = opt.lyap_skip; end
if isfield(opt, 'seed') && ~isempty(opt.seed)
    randn('state', opt.seed); rand('state', opt.seed);
end

N  = C.N; NE = P.NE;
JE = C.JE; JI = C.JI;
isE = C.isE;
sE = ones(N, 1); if isfield(opt, 'scaleE') && ~isempty(opt.scaleE), sE = opt.scaleE(:); end
gA = P.gA * sE;  gN = P.gN * sE;  gG = P.gG;
[kc, EAe, EGe] = syn_mode_coeffs(P);   % conductance- vs current-based AMPA/GABA

dt = P.dt; nsub = round(P.bin_ms / dt);
km = dt / P.tau_m;
eA = exp(-dt / P.tau_A); eG = exp(-dt / P.tau_G);
eN = exp(-dt ./ P.tau_N(:));   % tau_N may be a scalar or an NE x 1 vector (per presynaptic unit)
noise_amp = P.sig_eta * sqrt(dt / P.tau_m);

% initial state
if isfield(opt, 'x0') && ~isempty(opt.x0)
    V = opt.x0.V; a = opt.x0.a; n = opt.x0.n; b = opt.x0.b;
else
    V0 = P.E_L + 15;  if isfield(opt, 'V_init'), V0 = opt.V_init; end
    V = V0 + 2 * randn(N, 1);
    nu = rate_fn(V, P);
    a = nu(isE); n = nu(isE); b = nu(~isE);
end

useN = any(gN > 0);   % without NMDA the n variables are decoupled and excluded
if dolyap
    dV = randn(N, 1); da = randn(NE, 1); dn = randn(NE, 1) * useN; db = randn(N - NE, 1);
    nrm = sqrt(sum(dV.^2) + sum(da.^2) + sum(dn.^2) + sum(db.^2));
    dV = dV / nrm; da = da / nrm; dn = dn / nrm; db = db / nrm;
    logsum = 0; nacc = 0;
end

R = zeros(N, T_bins);
sumE = 0; sumI = 0;
Ibase = P.I_bg + P.f_in * P.I_ff;
for k = 1:T_bins
    Iin = Ibase + P.f_in * P.sig_u * C.w * u(k);
    for s = 1:nsub
        [nu, dnu] = rate_fn(V, P);
        [B, dB]   = mg_block(V, P);
        A  = JE * a;  Nn = JE * n;  G = JI * b;
        % exponential-Euler step for V (exact for conductances frozen over dt)
        Gt   = 1 + kc * (gA .* A + gG * G) + gN .* Nn .* B;        % total conductance
        U    = P.E_L + gA .* A * EAe + gN .* Nn .* B * P.E_E + gG * G * EGe + Iin;
        Vinf = U ./ Gt;
        Ee   = exp(-km * Gt);
        if dolyap
            % exact derivative of the discrete update (V, a, n, b) -> (V', a', n', b')
            dA = JE * da;  dNn = JE * dn;  dGg = JI * db;  dBv = dB .* dV;
            dNM  = gN .* (B .* dNn + Nn .* dBv);
            dGt  = kc * (gA .* dA + gG * dGg) + dNM;
            dU   = gA .* dA * EAe + dNM * P.E_E + gG * dGg * EGe;
            dVinf = (dU - Vinf .* dGt) ./ Gt;
            dVnew = (1 - Ee) .* dVinf + Ee .* dV - km * Ee .* (V - Vinf) .* dGt;
            dnuV = dnu .* dV;
            da = eA * da + (1 - eA) * dnuV(isE);
            dn = (eN .* dn + (1 - eN) .* dnuV(isE)) * useN;
            db = eG * db + (1 - eG) * dnuV(~isE);
            dV = dVnew;
        end
        V = Vinf + (V - Vinf) .* Ee;
        if noise, V = V + noise_amp * randn(N, 1); end
        a = nu(isE)  + (a - nu(isE))  * eA;
        n = nu(isE)  + (n - nu(isE))  .* eN;
        b = nu(~isE) + (b - nu(~isE)) * eG;
    end
    if any(~isfinite(V))
        error('simulate_net:unstable', 'Integration became unstable (reduce P.dt).');
    end
    nu = rate_fn(V, P);
    R(:, k) = nu;
    sumE = sumE + mean(nu(isE)); sumI = sumI + mean(nu(~isE));
    if dolyap
        nrm = sqrt(sum(dV.^2) + sum(da.^2) + sum(dn.^2) + sum(db.^2));
        if k > skip
            logsum = logsum + log(nrm); nacc = nacc + 1;
        end
        dV = dV / nrm; da = da / nrm; dn = dn / nrm; db = db / nrm;
    end
end

out.R = R;
out.x = struct('V', V, 'a', a, 'n', n, 'b', b);
out.rateE = P.r_max * sumE / T_bins;
out.rateI = P.r_max * sumI / T_bins;
if dolyap && nacc > 0
    out.lam = logsum / (nacc * P.bin_ms * 1e-3);   % 1/s
else
    out.lam = NaN;
end
end
