function [kc, EAe, EGe] = syn_mode_coeffs(P)
%SYN_MODE_COEFFS  Coefficients for conductance- or current-based AMPA/GABA.
%   P.syn_mode = 'conductance' (default): I_X = g_X s_X (E_X - V)
%   P.syn_mode = 'current'    : I_X = g_X s_X (E_X - V_ref), a fixed
%                               driving force evaluated at P.V_ref.
%   NMDA always keeps its full voltage dependence B(V)(E_E - V).
%   kc multiplies the AMPA/GABA conductances in the total conductance;
%   EAe, EGe are the effective reversal terms used in the update.
if ~isfield(P, 'syn_mode') || strcmp(P.syn_mode, 'conductance')
    kc = 1; EAe = P.E_E; EGe = P.E_I;
else
    Vr = -55; if isfield(P, 'V_ref'), Vr = P.V_ref; end
    kc = 0; EAe = P.E_E - Vr; EGe = P.E_I - Vr;
end
end
