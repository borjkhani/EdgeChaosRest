function C = build_connectivity(P, seed)
%BUILD_CONNECTIVITY  Sparse random E/I connectivity obeying Dale's law.
%   Each unit receives exactly KE excitatory and KI inhibitory inputs
%   (no self-connections). Weights are lognormal with mean 1/K and
%   coefficient of variation P.cv_w (before capping at P.w_max times the
%   mean, default 5), so that every row of JE and JI sums
%   to 1 on average. Feedforward input weights w are i.i.d. N(0,1).
%
%   Returns C.JE (N x NE), C.JI (N x NI), C.w (N x 1), C.N, C.isE.
if nargin > 1 && ~isempty(seed)
    rand('state', seed); randn('state', seed);   % works in MATLAB and Octave
end
N  = P.NE + P.NI;
s2 = log(1 + P.cv_w^2);
mu = -s2 / 2;
wmax = 5;                      % weights are capped at 5x the mean (heavy tail)
if isfield(P, 'w_max'), wmax = P.w_max; end
JE = zeros(N, P.NE);
JI = zeros(N, P.NI);
for i = 1:N
    candE = 1:P.NE;  candE(candE == i) = [];
    candE = candE(randperm(numel(candE)));
    idxE  = candE(1:P.KE);
    JE(i, idxE) = min(exp(mu + sqrt(s2) * randn(1, P.KE)), wmax) / P.KE;
    candI = 1:P.NI;  candI(candI == i - P.NE) = [];
    candI = candI(randperm(numel(candI)));
    idxI  = candI(1:P.KI);
    JI(i, idxI) = min(exp(mu + sqrt(s2) * randn(1, P.KI)), wmax) / P.KI;
end
C.JE  = JE;
C.JI  = JI;
C.w   = randn(N, 1);
C.N   = N;
C.isE = [true(P.NE, 1); false(P.NI, 1)];
end
