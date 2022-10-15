function [gammaBar_dB, Pb] = OP_analit(alpha, mu, ms, bounds, N, z, gamma_th)
% TODO: doc func
% N - number of points

%lambda = 1;

% generate independent variable vector
L = bounds(1);
U = bounds(2);
gammaBar_dB = linspace(L, U, N);
%gammaBar = gpuArray.linspace(L, U, N);

gammaBar = db2pow(gammaBar_dB);

uPsi = (mu/(ms-1)) * ((z*sqrt(gamma_th)) ./ (sqrt(gammaBar*(z^2+2)))) .^ alpha;

% precomputations
preGammaCoef = z^2/(alpha*gamma(mu)*gamma(ms));

% m = 2; n = 2; p = 3; q = 3;
an = [(1-ms), 1];                       An = [1, 1];   
ap = [((z^2)/alpha)+1];                 Ap = [1];
bm = [mu,((z^2)/alpha)];                Bm = [1,1];
bq = [0];                               Bq = [1];


%fox = zeros(1, N, 'gpuArray');
fox = zeros(1, N);
for i = 1:N
    fox(i) = real(HFox(an, An, ap, Ap, bm, Bm ,bq, Bq, uPsi(i)));
    %fox(i) = real(gpuArray(HFox(an, An, ap, Ap, bm, Bm ,bq, Bq, uPsi(i))));
end

% compute OP
Pb = preGammaCoef .* fox;

% DEBUG
% semilogy(gammaBar_dB, Pb)
end