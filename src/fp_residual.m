function [F, JF] = fp_residual(P, C, V, scaleE)
%FP_RESIDUAL  Fixed-point equations F(V) = 0 of the network with constant
%   input and their Jacobian dF/dV = -diag(G) + M(0) Phi' (see newton_fp.m).
if nargin < 4 || isempty(scaleE), scaleE = ones(C.N, 1); end
isE = C.isE; gA = P.gA * scaleE; gN = P.gN * scaleE; gG = P.gG;
[kc, EAe, EGe] = syn_mode_coeffs(P);
I = P.I_bg + P.f_in * P.I_ff;
[nu, dnu] = rate_fn(V, P); [B, dB] = mg_block(V, P);
A = C.JE * nu(isE); G = C.JI * nu(~isE); DN = P.E_E - V;
if kc, DA = P.E_E - V; DI = P.E_I - V; else, DA = EAe; DI = EGe; end
F = -(V - P.E_L) + gA .* A .* DA + gN .* B .* A .* DN + gG * G .* DI + I;
if nargout > 1
    Gloc = 1 + kc * (gA .* A + gG * G) + gN .* A .* (B - dB .* DN);
    rowE = gA .* DA + gN .* B .* DN; rowI = gG * DI;
    JF = -diag(Gloc) + [diag(rowE) * C.JE * diag(dnu(isE)), diag(rowI) * C.JI * diag(dnu(~isE))];
end
end
