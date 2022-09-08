function [gammaBar_dB, Pb] = SEP_analit(theta, alpha, betaVar, mu, ms, phi, bounds, N, hl, z)
% TODO: doc func
% N - number of points

%lambda = 1;

% compute Lambda0
A_0 = sqrt( gamma(3/betaVar) / gamma(1/betaVar) );
%A_0 = 0.8;

% generate independent variable vector
L = bounds(1);
U = bounds(2);
gammaBar = linspace(L, U, N);
%gammaBar = gpuArray.linspace(L, U, N);

gammaBar_dB = pow2db(gammaBar);

%uPsi = ( mu./(ms-1)*(sqrt(2))./ sqrt(phi).*A_0 .* gammaBar.^(1/2) .* hl );
uPsi = (mu/(ms-1))*((sqrt(2)./(sqrt(phi)*A_0 * gammaBar.^(1/2)*hl)).^alpha);

% precomputations
%preGammaCoef = theta / (betaVar * gamma(1/betaVar) * gamma(mu) * gamfma(ms) );
preGammaCoef = (theta / (2*sqrt(pi))) * (z^2 / (alpha * gamma(mu) * gamma(ms))); 

% m = 2; n = 3; p = 4; q = 3;
an = [(1-ms), 1, 1/2];                  An = [1, 1, (alpha/2)];   
ap = [((z^2)/alpha)+1];                 Ap = [1];
bm = [mu,((z^2)/alpha)];                Bm = [1,1];
bq = [0];                               Bq = [1];


%fox = zeros(1, N, 'gpuArray');
fox = zeros(1, N);
for i = 1:N
    fox(i) = real(HFox(an, An, ap, Ap, bm, Bm ,bq, Bq, uPsi(i)));
    %fox(i) = real(gpuArray(HFox(an, An, ap, Ap, bm, Bm ,bq, Bq, uPsi(i))));
end

% compute SEP
Pb = preGammaCoef .* fox;

% DEBUG
% semilogy(gammaBar_dB, Pb)
end
