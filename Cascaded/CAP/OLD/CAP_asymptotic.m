function [gammaBar_dB, P] = CAP_asymptotic(N, alpha, mu, ms, bounds, points, z)
% Function to implement the asymptotic expression 
% for the alpha-F fading and AWGGN distributions.
%
% INPUTS:
% alpha - positive power parameter
% beta - 
% mu - number of multipath clusters
% ms - shadowing parameter
% points - number of points
%
% OUTPUTS:
% P - points that form the asymptote

% generate independent variable vector
L = bounds(1);
U = bounds(2);
gammaBar_dB = linspace(L, U, points);
%gammaBar = gpuArray.linspace(L, U, points);

gammaBar = db2pow(gammaBar_dB);

Xi = ones(1, length(gammaBar));
S = 0;
for j = 1:N
    Psi = (mu/(ms-1)) ^ (1/alpha);
    Xi = Xi .* (Psi .* (z ./ (sqrt(gammaBar .* (z^2+2)))));

    frac = (z^2) / alpha;
    S = S + (2/alpha)*(psi(mu) + psi(frac) - psi(ms) - psi(frac+1));
end

P = (2*log2(1 ./ Xi)) + (log2(exp(1))*S);

end