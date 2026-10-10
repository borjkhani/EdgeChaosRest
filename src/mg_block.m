function [B, dB] = mg_block(V, P)
%MG_BLOCK  Voltage-dependent Mg2+ unblock of NMDA receptors (Jahr & Stevens 1990).
%   B(V)  = 1 ./ (1 + [Mg]/gamma * exp(-delta*V))
%   dB/dV = delta * B .* (1 - B)
B  = 1 ./ (1 + (P.Mg / P.mg_gamma) .* exp(-P.mg_delta .* V));
dB = P.mg_delta .* B .* (1 - B);
end
