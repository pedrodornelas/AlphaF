function [gammaBar_dB, Pb] = CAP_analit(N, alpha, mu, ms, bounds, points, z)
% TODO: doc func
% points - number of points

%lambda = 1;

% generate independent variable vector
L = bounds(1);
U = bounds(2);
gammaBar_dB = linspace(L, U, points);
%gammaBar = gpuArray.linspace(L, U, points);

gammaBar = db2pow(gammaBar_dB);

Xi = ones(1, length(gammaBar));
preGammaCoef = 1;
for j = 1:N
    Psi = (mu/(ms-1)) ^ (1/alpha);
    Xi = Xi .* (Psi .* (z ./ (sqrt(gammaBar .* (z^2+2)))));

    % precomputations
    preGammaCoef = preGammaCoef * (z^2/(alpha*gamma(mu)*gamma(ms)));
end

preGammaCoef = preGammaCoef / (2*log(2));

onesN = ones(1, N);

% m = 2N; n = N+2; p = 2N+2; q = 2N+1;
an = [0, (1-ms)*onesN];                      An = [1/2, (1/alpha)*onesN];   
ap = [((z^2)/alpha)*onesN + 1 , 1];          Ap = [(1/alpha)*onesN, 1/2];
bm = [mu*onesN, ((z^2)/alpha)*onesN, 0, 0];  Bm = [(1/alpha)*onesN, (1/alpha)*onesN, 1/2, 1/2];
bq = [];                                     Bq = [];


%fox = zeros(1, points, 'gpuArray');
fox = zeros(1, points);
for i = 1:points
    fox(i) = real(HFox(an, An, ap, Ap, bm, Bm ,bq, Bq, Xi(i)));
    %fox(i) = real(gpuArray(HFox(an, An, ap, Ap, bm, Bm ,bq, Bq, Xi(i))));
end

% compute CAP
Pb = preGammaCoef .* fox;

% DEBUG
% semilogy(gammaBar_dB, Pb)
end