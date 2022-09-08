function [gammaBar_dB, P] = OP_asymptotic(alpha, betaVar, mu, ms, bounds, N, hl, z, gamma_th)
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

% compute Lambda0
A_0 = sqrt( gamma(3/betaVar) / gamma(1/betaVar) );
%A_0 = 0.8;

Psi = mu/(ms-1);

% precomputations
%preBeta = theta * gamma( (1+alpha*mu)/betaVar) / 2 / mu / beta(mu, ms) / gamma(1/betaVar);
%prePow = ( (mu*lambda^(alpha/2)) / ( (A_0*sqrt(phi))^alpha * (ms-1)) )^mu;
preBeta = 0;
%prePow = (sqrt(phi)*hl*A_0*gammaBar.^(1/2)) / sqrt(2);
prePow = (hl*A_0*sqrt(gammaBar))/sqrt(gamma_th);
aux = 0;

% asymptotic OP
if mu < (z^2/alpha)
%    preBeta = (theta * z^2 * gamma((mu*alpha+1)/2) * Psi^mu) / ...
%              (2*sqrt(pi)*mu*(z^2 - mu*alpha)*beta(mu, ms));

    preBeta = (z.^2 * Psi.^mu)/((z.^2*mu - alpha*mu.^2) * (beta(mu,ms)));
    aux = -mu*alpha;
else
    preBeta = (gamma(mu-z.^2/alpha) * gamma(ms+z.^2/alpha)*Psi.^(z.^2/alpha))/(gamma(mu)*gamma(ms));
    aux = -z^2;
end
P = preBeta * prePow.^(aux);

% debug
% semilogy(gammaBar_dB, P)

end
