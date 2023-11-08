function [CAP] = CAP_asymptotic(N, params, gammaBar)
% params is a matrix of params cascaded system
% params = [Channel 1: alpha, mu, ms, z;
%           Channel 2: alpha, mu, ms, z;
%           ...
%           Channel n: alpha, mu, ms, z]

% alpha: non-linearity
% mu: number of multipaths
% ms: shadowing
% z: pointing error
% gamma_th: OP treshold
% gammaBar: SNR matriz per channel
% gammaBar = [          , channel1, channel2, ..., channeln;
%              gammaBar1,     y1.1,     y1.2, ...,     y1.n;
%                    ...,      ...,      ..., ...,      ...;
%             gammaBar15,    y15.1,    y15.2, ...,    y15.n;]

channels = N;
points = length(gammaBar);

Xi = 1;
S = 0;
for c = 1:N

    alpha = params(c,1);
    mu = params(c,2);
    ms = params(c,3);
    z = params(c,4);

    Psi = (mu/(ms-1)) ^ (1/alpha);
    Xi = Xi .* ((Psi .* z) ./ sqrt(gammaBar(:, c) .* (z^2+2)));

    frac = (z^2) / alpha;
    S = S + (2/alpha)*(psi(mu) + psi(frac) - psi(ms) - psi(frac+1));
end

CAP = (2*log2(1 ./ Xi)) + (log2(exp(1))*S);

end