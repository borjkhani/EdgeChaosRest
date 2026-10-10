function P = condition_params(k)
%CONDITION_PARAMS  The four synaptic conditions of the visual-cortex experiments.
%   1 AMPA only (fN = 0); 2 fN = 0.5, [Mg2+] = 1 mM; 3 fN = 0.5, [Mg2+] = 2 mM;
%   4 fN = 0.5 with 50% NMDA conductance block (memantine-like).
%   All at total excitatory conductance gE_tot = 1.4 (healthy, ordered).
gE = 1.4;
switch k
    case 1, P = set_composition(base_params(), 0,   gE, -55);
    case 2, P = set_composition(base_params(), 0.5, gE, -55);
    case 3, P = set_composition(base_params(), 0.5, gE, -55); P.Mg = 2;
    case 4, P = set_composition(base_params(), 0.5, gE, -55); P.gN = 0.5 * P.gN;
end
P.cond = k;
end
