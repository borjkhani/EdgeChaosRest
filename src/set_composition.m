function P = set_composition(P, fN, gE_tot, V_ref)
%SET_COMPOSITION  Split a fixed excitatory conductance between AMPA and NMDA.
%   fN is the NMDA fraction of the excitatory steady-state conductance
%   evaluated at the reference voltage V_ref (mV):
%       gA = (1 - fN) * gE_tot
%       gN * B(V_ref) = fN * gE_tot      =>  gN = fN * gE_tot / B(V_ref)
%   so that, at V = V_ref, the total excitatory synaptic conductance (and
%   hence the mean excitatory drive) does not depend on fN. Changing fN
%   therefore changes only (i) the voltage dependence of the excitatory
%   conductance (Mg2+ block) and (ii) the fraction carried by the slow
%   (tau_N) rather than the fast (tau_A) synapse.
if nargin < 4 || isempty(V_ref), V_ref = -55; end
Bref = mg_block(V_ref, P);
P.gA = (1 - fN) * gE_tot;
P.gN = fN * gE_tot / Bref;
P.fN = fN; P.gE_tot = gE_tot; P.V_ref = V_ref;
end
