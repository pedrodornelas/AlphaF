function [gammaBar_dB, Pb] = Capacity_analit( alpha, A_0, mu, ms, bounds, N, hl, z)
% TODO: doc func
% N - number of points

%lambda = 1;

% generate independent variable vector
L = bounds(1);
U = bounds(2);
gammaBar = linspace(L, U, N);
%gammaBar = gpuArray.linspace(L, U, N);

gammaBar_dB = pow2db(gammaBar);

uPsi = (mu/(ms-1))./((A_0.*sqrt(gammaBar)*hl).^alpha);

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
