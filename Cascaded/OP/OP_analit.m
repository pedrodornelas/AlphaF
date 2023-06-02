function [gammaBar_dB, Pb] = OP_analit(N, alpha, mu, ms, bounds, points, z, gamma_th)
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
    Psi = (mu(j)/(ms-1)) ^ (1/alpha);
    Xi = Xi .* (Psi .* (z(j) ./ (sqrt(gammaBar .* (z(j)^2+2)))));

    % precomputations
    preH = preH * (z(j)^2/((alpha)*gamma(mu(j))*gamma(ms)));
end

Xi = Xi * sqrt(gamma_th);
onesAn = ones(1, N);

% m = 2N; n = N+1; p = 2N+1; q = 2N+1;
an = [1, (1-ms)*onesAn];                 An = [1, (1/alpha)*onesAn];   
ap = [(((z.^2)./alpha) + 1)];            Ap = [(1/alpha)*onesAn];
bm = [mu, ((z.^2)./alpha)];              Bm = [(1/alpha)*onesAn, (1/alpha)*onesAn];
bq = [0];                                Bq = [1];


%fox = zeros(1, points, 'gpuArray');
fox = zeros(1, points);
for i = 1:points
    fox(i) = real(HFox(an, An, ap, Ap, bm, Bm ,bq, Bq, Xi(i)));
    %fox(i) = real(gpuArray(HFox(an, An, ap, Ap, bm, Bm ,bq, Bq, Xi(i))));
end

% compute OP
Pb = preH .* fox;

% DEBUG
% semilogy(gammaBar_dB, Pb)
end