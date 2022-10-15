function [gammaBar_dB, P] = OP_asymptotic(alpha, mu, ms, bounds, N, z, gamma_th)
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
gammaBar_dB = linspace(L, U, N);
%gammaBar = gpuArray.linspace(L, U, N);

gammaBar = db2pow(gammaBar_dB);

Psi = mu/(ms-1);

% precomputations
preBeta = 0;
prePow = (sqrt(z^2+2)*sqrt(gammaBar))/(z*sqrt(gamma_th));
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