function [gammaBar_dB, Pb] = Capacity_analit(alpha, mu, ms, bounds, N, z)
% TODO: doc func
% N - number of points

% generate independent variable vector
L = bounds(1);
U = bounds(2);
gammaBar_dB = linspace(L, U, N);
%gammaBar = gpuArray.linspace(L, U, N);

gammaBar = db2pow(gammaBar_dB);

uPsi = (mu/(ms-1)) * ((z ./ (sqrt(z^2+2)*sqrt(gammaBar))).^alpha);

% precomputations

preGammaCoef = z^2/(2*log(2)*gamma(mu)*gamma(ms));

% m = 4; n = 2; p = 4; q = 4;
an = [(1-ms), 0];                       An = [1, (alpha/2)];   
ap = [1,((z^2)/alpha +1)];             Ap = [(alpha/2),1];
bm = [mu,((z^2)/alpha), 0, 0];          Bm = [1,1,alpha/2,alpha/2];
bq = [];                               Bq = [];


%fox = zeros(1, N, 'gpuArray');
fox = zeros(1, N);
for i = 1:N
    fox(i) = real(HFox(an, An, ap, Ap, bm, Bm ,bq, Bq, uPsi(i)));
    %fox(i) = real(gpuArray(HFox(an, An, ap, Ap, bm, Bm ,bq, Bq, uPsi(i))));
end

% compute BEP
Pb = preGammaCoef .* fox;

% DEBUG
% semilogy(gammaBar_dB, Pb)
end