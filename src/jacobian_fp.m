function [Jac, info] = jacobian_fp(P, C, x, scaleE)
%JACOBIAN_FP  Jacobian of the continuous-time network at state x.
%   State ordering: [V (N); a (NE); n (NE); b (NI)]. Time unit: ms.
%   info returns the ingredients used in the analytic gain decomposition:
%     info.Gloc   : local slope -dF_i/dV_i / 1 (total effective conductance)
%     info.Gnsc   : NMDA negative-slope-conductance contribution (<0 when
%                   the Mg2+ term dominates)
%     info.dnu    : phi'(V)
N = C.N; NE = P.NE; NI = P.NI; isE = C.isE;
if nargin < 4 || isempty(scaleE), scaleE = ones(N, 1); end
gA = P.gA * scaleE; gN = P.gN * scaleE; gG = P.gG;
V = x.V; a = x.a; n = x.n; b = x.b;
[~, dnu] = rate_fn(V, P);
[B, dB]  = mg_block(V, P);
A = C.JE * a; Nn = C.JE * n; G = C.JI * b;
[kc, EAe, EGe] = syn_mode_coeffs(P);
dfN = P.E_E - V;                                     % NMDA driving force (always voltage dependent)
if kc, dfE = P.E_E - V; dfI = P.E_I - V; else, dfE = EAe * ones(N, 1); dfI = EGe * ones(N, 1); end
Gnsc = -gN .* Nn .* dB .* dfN;                       % = gN*Nn*dB*(V-E_E) <= 0
Gloc = 1 + kc * (gA .* A + gG * G) + gN .* Nn .* B + Gnsc;  % total slope conductance
tm = P.tau_m;
JVV = -diag(Gloc) / tm;
JVa = diag(gA .* dfE) * C.JE / tm;
JVn = diag(gN .* B .* dfN) * C.JE / tm;
JVb = diag(gG * dfI) * C.JI / tm;
D = diag(dnu);
Jaa = -eye(NE) / P.tau_A; JaV = D(isE, :) / P.tau_A;
tN = P.tau_N(:) .* ones(NE, 1);                     % scalar or per-unit NMDA decay
Jnn = -diag(1 ./ tN); JnV = diag(1 ./ tN) * D(isE, :);
Jbb = -eye(NI) / P.tau_G; JbV = D(~isE, :) / P.tau_G;
Z = @(r, c) zeros(r, c);
Jac = [JVV, JVa, JVn, JVb; ...
       JaV, Jaa, Z(NE, NE), Z(NE, NI); ...
       JnV, Z(NE, NE), Jnn, Z(NE, NI); ...
       JbV, Z(NI, NE), Z(NI, NE), Jbb];
info.Gloc = Gloc; info.Gnsc = Gnsc; info.dnu = dnu; info.B = B;
info.A = A; info.Nn = Nn; info.G = G;
info.dfE = dfE; info.dfI = dfI; info.dfN = dfN; info.kc = kc;
end
