function S = gain_criterion(P, C, x, scaleE)
%GAIN_CRITERION  Static-gain criterion for the instability of a fixed point.
%
%   Linearising the network at a fixed point x and eliminating the synaptic
%   gating variables in the Laplace domain gives, for dV ~ exp(s t),
%
%       (tau_m s + G_i) dV_i = sum_j M_ij(s) phi'_j dV_j ,
%
%   with the local slope conductance (dimensionless, units of g_L)
%       G_i = 1 + gA A_i + gN N_i [B(V_i) + B'(V_i)(V_i - E_E)] + gG Gi_i
%   (the B' term is the NMDA negative slope conductance, NSC <= 0 for V<E_E)
%   and the synaptic coupling
%       M_ij(s) = (E_E - V_i)[gA H_A(s) + gN B(V_i) H_N(s)] J^E_ij  (j in E)
%       M_ij(s) = (E_I - V_i) gG H_G(s) J^I_ij                      (j in I)
%   with H_X(s) = 1/(1 + tau_X s). At s = 0 all filters equal 1, so the
%   Jacobian is singular exactly when the static gain matrix
%       K = diag(1./G) * M(0) * diag(phi')
%   has an eigenvalue equal to 1. K contains no time constant.
%
%   Outputs:
%     S.rhoK        : largest real part of eig(K)   (static gain)
%     S.rhoK_noNSC  : same with the NSC removed from G (ablation)
%     S.rhoK_noMg   : same with NMDA treated as voltage-independent,
%                     B(V) = B(V_ref), no NSC (ablation of the Mg2+ gate)
%     S.maxReJ      : largest real part of eig(Jacobian) (1/s)
%     S.leadComplex : true if the leading Jacobian eigenvalue is complex
%     S.G, S.NSC, S.V, S.dnu : population means at the fixed point
if nargin < 4 || isempty(scaleE), scaleE = ones(C.N, 1); end
[Jac, info] = jacobian_fp(P, C, x, scaleE);
gA = P.gA * scaleE; gN = P.gN * scaleE;
V = x.V; isE = C.isE;
B = mg_block(V, P);
dnu = info.dnu;
G     = info.Gloc;
G_noN = G - info.Gnsc;
DA = info.dfE; DN = info.dfN; DI = info.dfI;
rowE = gA .* DA + gN .* B .* DN;
S.rhoK       = static_gain(G,     rowE, P.gG * DI, C, dnu, isE);
S.rhoK_noNSC = static_gain(G_noN, rowE, P.gG * DI, C, dnu, isE);
Vr = -55; if isfield(P, 'V_ref'), Vr = P.V_ref; end
Bref = mg_block(Vr, P);
G_noMg = 1 + info.kc * (gA .* info.A + P.gG * info.G) + gN .* info.Nn * Bref;
S.rhoK_noMg  = static_gain(G_noMg, gA .* DA + gN * Bref .* DN, P.gG * DI, C, dnu, isE);
ev = eig(Jac);
[~, im] = max(real(ev));
S.maxReJ = 1000 * real(ev(im));
S.leadComplex = abs(imag(ev(im))) > 1e-9;
S.G = mean(G); S.NSC = mean(info.Gnsc); S.V = mean(V); S.dnu = mean(dnu);
end

function r = static_gain(G, rowE, rowI, C, dnu, isE)
K = [diag(rowE ./ G) * C.JE * diag(dnu(isE)), diag(rowI ./ G) * C.JI * diag(dnu(~isE))];
r = max(real(eig(K)));
end
