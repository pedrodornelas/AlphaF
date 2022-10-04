function [gammaBar_dB, P] = SEP_asymptotic(theta,alpha, mu, ms, phi, bounds, N, z)
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

lambda = 1;

% generate independent variable vector
L = bounds(1);
U = bounds(2);
%gammaBar = gpuArray.linspace(L, U, N);
gammaBar = linspace(L, U, N);

gammaBar_dB = pow2db(gammaBar);

Psi = mu/(ms-1);

% precomputations
%preBeta = theta * gamma( (1+alpha*mu)/betaVar) / 2 / mu / beta(mu, ms) / gamma(1/betaVar);
%prePow = ( (mu*lambda^(alpha/2)) / ( (A_0*sqrt(phi))^alpha * (ms-1)) )^mu;
preBeta = 0;
prePow = (sqrt(phi*(z^2+2)) * gammaBar.^(1/2)) / (z*sqrt(2));
aux = 0;

% asymptotic SEP
if mu < (z^2/alpha)
    preBeta = (theta * z^2 * gamma((mu*alpha+1)/2) * Psi^mu) / ...
              (2*sqrt(pi)*mu*(z^2 - mu*alpha)*beta(mu, ms));
    aux = -mu*alpha;
else
    preBeta = (theta * gamma(mu-((z^2)/alpha)) * gamma(ms+z^2/alpha) * gamma((1+z^2)/2) * Psi.^(z^2/alpha)) / ...
              (2*sqrt(pi)*gamma(mu)*gamma(ms));
    aux = -z^2;
end
P = preBeta * prePow.^(aux);

% debug
% semilogy(gammaBar_dB, P)

end
