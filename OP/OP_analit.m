function [gammaBar_dB, Pb] = OP_analit( alpha, betaVar, mu, ms, bounds, N, hl, z, gamma_th)
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

%uPsi = (mu/(ms-1))*((sqrt(2)./(sqrt(phi)*A_0 * gammaBar.^(1/2)*hl)).^alpha);
uPsi = (mu/(ms-1)) * ((sqrt(gamma_th))./(sqrt(gammaBar)*hl*A_0)).^alpha;
uPsi = (mu/(ms-1)) * ((sqrt(gamma_th)) ./ (sqrt(gammaBar) * hl*A_0)) .^ alpha;

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
