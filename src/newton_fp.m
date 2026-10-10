function [x, ok, info] = newton_fp(P, C, V0, scaleE, tol, maxit)
%NEWTON_FP  Fixed point of the network by Newton's method (stable or not).
%   At a fixed point every gating variable equals phi(V), so the fixed
%   points are the roots of F(V) = 0 with
%     F_i = -(V_i - E_L) + theta gA A_i DA_i + theta gN B(V_i) A_i (E_E - V_i)
%           + gG Gi_i DI_i + I_i,   A = J^E phi(V_E),  Gi = J^I phi(V_I),
%   (DA, DI: AMPA/GABA driving forces; voltage dependent for conductance-
%   based synapses). dF/dV = -diag(G) + M(0) diag(phi'), whose singularity is
%   the static-gain criterion of gain_criterion.m. No time constant enters.
%   Input: constant drive I_bg + f_in*I_ff (u = 0, no noise).
%   Returns the full state x (V, a = n = phi(V_E), b = phi(V_I)), a success
%   flag and info.normF, info.iter.
if nargin < 4 || isempty(scaleE), scaleE = ones(C.N, 1); end
if nargin < 5 || isempty(tol), tol = 1e-9; end
if nargin < 6 || isempty(maxit), maxit = 60; end
isE = C.isE; gA = P.gA * scaleE; gN = P.gN * scaleE; gG = P.gG;
[kc, EAe, EGe] = syn_mode_coeffs(P);
I = P.I_bg + P.f_in * P.I_ff;
V = V0(:); ok = false;
for it = 1:maxit
    [nu, dnu] = rate_fn(V, P); [B, dB] = mg_block(V, P);
    A = C.JE * nu(isE); G = C.JI * nu(~isE);
    DN = P.E_E - V;
    if kc, DA = P.E_E - V; DI = P.E_I - V; else, DA = EAe; DI = EGe; end
    F = -(V - P.E_L) + gA .* A .* DA + gN .* B .* A .* DN + gG * G .* DI + I;
    nrm = max(abs(F));
    if nrm < tol, ok = true; break; end
    Gloc = 1 + kc * (gA .* A + gG * G) + gN .* A .* (B - dB .* DN);
    rowE = gA .* DA + gN .* B .* DN; rowI = gG * DI;
    JF = -diag(Gloc) + [diag(rowE) * C.JE * diag(dnu(isE)), diag(rowI) * C.JI * diag(dnu(~isE))];
    dV = -JF \ F;
    if any(~isfinite(dV)), break; end
    step = 1; % damping: keep the voltage update bounded
    mx = max(abs(dV)); if mx > 5, step = 5 / mx; end
    V = V + step * dV;
end
nu = rate_fn(V, P);
x = struct('V', V, 'a', nu(isE), 'n', nu(isE), 'b', nu(~isE));
info = struct('normF', nrm, 'iter', it);
end
