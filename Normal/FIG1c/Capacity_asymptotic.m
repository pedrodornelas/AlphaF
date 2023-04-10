function [gammaBar_dB, P] = Capacity_asymptotic(alpha, A_0, mu, ms, bounds, N, hl, z)
% Function to implement the asymptotic expression 
% for the alpha-F fading and AWGGN distributions.
%
% INPUTS:
% alpha - positive power parameter
% beta - 
% mu - number of multipath clusters
% ms - shadowing parameter
% phi -
% N - number of points
%
% OUTPUTS:
% P - points that form the asymptote

% generate independent variable vector
L = bounds(1);
U = bounds(2);
%gammaBar = gpuArray.linspace(L, U, N);
gammaBar = linspace(L, U, N);

gammaBar_dB = pow2db(gammaBar);

P = log2(gammaBar)+2*log2(exp(1))* ...
    (psi(mu)/alpha - (psi(ms)/(alpha)) - 1/(z^2) + log(((ms-1)/mu)^(1/alpha)*A_0*hl));

% debug
% semilogy(gammaBar_dB, P)

end
