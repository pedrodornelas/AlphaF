function [gammaBar_dB, P] = OP_asymptotic(N, alpha, mu, ms, bounds, points, z, gamma_th)
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

lambda = 1;

% generate independent variable vector
L = bounds(1);
U = bounds(2);
gammaBar_dB = linspace(L, U, points);
%gammaBar = gpuArray.linspace(L, U, points);

gammaBar = db2pow(gammaBar_dB);

Xi = ones(1, length(gammaBar));
preGammaCoef = 1;
for j = 1:N
    Psi = (mu(j)/(ms-1)) ^ (1/alpha);
    Xi = Xi .* (Psi .* ((z(j).*sqrt(gamma_th)) ./ (sqrt(gammaBar .* (z(j)^2+2)))));

    % precomputations
    preGammaCoef = preGammaCoef * (z(j)^2/((alpha)*gamma(mu(j))*gamma(ms)));
end

onesAn = ones(1, N);

% m = 2N; n = N+1; p = 2N+1; q = 2N+1;
an = [1, (1-ms)*onesAn];                 An = [1, (1/alpha)*onesAn];   
ap = [(((z.^2)./alpha) + 1)];            Ap = [(1/alpha)*onesAn];
bm = [mu, ((z.^2)./alpha)];              Bm = [(1/alpha)*onesAn, (1/alpha)*onesAn];
bq = [0];                                Bq = [1];


div = bm./Bm;
[U, idx] = min(div);

mins = find(div == U);

Bc = 1;
Bc = Bm(idx);

% gamma 1
gamma1 = 1;
for j = 1:length(bm)
    % j not is member of mins
    if ismember(j,mins) == 0
        gamma1 = gamma1*gamma(bm(j) - U*Bm(j));
    else
        % Bc = Bc*Bm(j);
    end
end

gamma2 = 1;
for j = 1:length(an)
    gamma2 = gamma2*gamma(1 - an(j) + U*An(j));
end

gamma3 = 1;
for j = 1:length(bq)
    gamma3 = gamma3*gamma(1 - bq(j) + U*Bq(j));
end

gamma4 = 1;
for j = 1:length(ap)
    gamma4 = gamma4*gamma(ap(j) - U*Ap(j));
end

phiU = (gamma1*gamma2) / (gamma3*gamma4);

P = (preGammaCoef * phiU * (Xi .^ U)) ./ Bc;


% precomputations
% preBeta = 0;
% prePow = (sqrt(z^2+2)*sqrt(gammaBar))/(z*sqrt(gamma_th));
% aux = 0;

% % asymptotic OP
% if mu < (z^2/alpha)
% %    preBeta = (theta * z^2 * gamma((mu*alpha+1)/2) * Psi^mu) / ...
% %              (2*sqrt(pi)*mu*(z^2 - mu*alpha)*beta(mu, ms));

%     preBeta = (z.^2 * Psi.^mu)/((z.^2*mu - alpha*mu.^2) * (beta(mu,ms)));
%     aux = -mu*alpha;
% else
%     preBeta = (gamma(mu-z.^2/alpha) * gamma(ms+z.^2/alpha)*Psi.^(z.^2/alpha))/(gamma(mu)*gamma(ms));
%     aux = -z^2;
% end
% P = preBeta * prePow.^(aux);

% debug
% semilogy(gammaBar_dB, P)

end