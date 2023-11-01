function [gammaBar_dB, Pb] = BEP_analit(N, alpha, mu, ms, bounds, points, z, rho)
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
preH = 1;
for j = 1:N
    Psi = (mu/(ms-1)) ^ (1/alpha(j));
    Xi = Xi .* (Psi .* (z(j) ./ (sqrt(gammaBar .* (z(j)^2+2)))));

    % precomputations
    preH = preH * (z(j)^2/(alpha(j)*gamma(mu)*gamma(ms)));
end

Xi = Xi * (1/sqrt(rho));
preH = preH / (4*sqrt(pi));

onesN = ones(1, N);

% m = 2N; n = N+2; p = 2N+2; q = 2N+1;
an = [(1-ms)*onesN, 1, (1/2)];           An = [(1./alpha), (1/2), (1/2)];   
ap = [(((z.^2)./alpha) + 1)];            Ap = [(1./alpha)];
bm = [mu*onesN, ((z.^2)./alpha)];        Bm = [(1./alpha), (1./alpha)];
bq = [0];                                Bq = [1/2];


%fox = zeros(1, points, 'gpuArray');
fox = zeros(1, points);
for i = 1:points
    fox(i) = real(HFox(an, An, ap, Ap, bm, Bm ,bq, Bq, Xi(i)));
    %fox(i) = real(gpuArray(HFox(an, An, ap, Ap, bm, Bm ,bq, Bq, Xi(i))));
end

% compute BEP
Pb = preH .* fox;

% DEBUG
% semilogy(gammaBar_dB, Pb)
end