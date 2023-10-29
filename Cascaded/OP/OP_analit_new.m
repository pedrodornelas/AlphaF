function [OP] = OP_analit_new(N, params, gamma_th, gammaBar)
% params is a matrix of params cascaded system
% params = [Channel 1: alpha, mu, ms, z;
%           Channel 2: alpha, mu, ms, z;
%           ...
%           Channel n: alpha, mu, ms, z]

% alpha: non-linearity
% mu: number of multipaths
% ms: shadowing
% z: pointing error
% gamma_th: OP treshold
% gammaBar: SNR vector

channels = N;
points = length(gammaBar);

Xi = 1;
preH = 1;
for c = 1:N

    alpha = params(c,1);
    mu = params(c,2);
    ms = params(c,3);
    z = params(c,4);

    Psi = (mu/(ms-1)) ^ (1/alpha);
    Xi = Xi * ((Psi * z) / sqrt(z^2+2));

    % precomputations
    preH = preH * (z^2/(alpha * gamma(mu) * gamma(ms)));
end

Xi = (Xi * sqrt(gamma_th)) ./ sqrt(gammaBar);

% m = 2N; n = N+1; p = 2N+1; q = 2N+1;
alphas = params(1:N, 1).';
mus = params(1:N, 2).';
mss = params(1:N, 3).';
zs = params(1:N, 4).';

an = [1, (1 - mss)];                 An = [1, (1./alphas)];   
ap = [(((zs.^2)./alphas) + 1)];      Ap = [(1./alphas)];
bm = [mus, ((zs.^2)./alphas)];       Bm = [(1./alphas), (1./alphas)];
bq = [0];                            Bq = [1];


%fox = zeros(1, points, 'gpuArray');
fox = zeros(1, points);
for i = 1:points
    fox(i) = real(HFox(an, An, ap, Ap, bm, Bm ,bq, Bq, Xi(i)));
    %fox(i) = real(gpuArray(HFox(an, An, ap, Ap, bm, Bm ,bq, Bq, Xi(i))));
end

% compute OP
OP = preH .* fox;

% DEBUG
% semilogy(gammaBar_dB, Pb)
end