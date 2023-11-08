function [CAP] = CAP_analit(N, params, gammaBar)
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
% gammaBar: SNR matriz per channel
% gammaBar = [          , channel1, channel2, ..., channeln;
%              gammaBar1,     y1.1,     y1.2, ...,     y1.n;
%                    ...,      ...,      ..., ...,      ...;
%             gammaBar15,    y15.1,    y15.2, ...,    y15.n;]

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
    Xi = Xi .* ((Psi .* z) ./ sqrt(gammaBar(:, c) .* (z^2+2)));

    % precomputations
    preH = preH * (z^2/(alpha * gamma(mu) * gamma(ms)));
end

preH = preH / (2*log(2));

alphas = params(1:N, 1).';
mus = params(1:N, 2).';
mss = params(1:N, 3).';
zs = params(1:N, 4).';

% m = 2N; n = N+2; p = 2N+2; q = 2N+1;
an = [0, (1-mss)];                   An = [1/2, (1./alphas)];   
ap = [(((zs.^2)./alphas) + 1), 1];   Ap = [(1./alphas), 1/2];
bm = [mus, ((z.^2)./alphas), 0, 0];  Bm = [(1./alphas), (1./alphas), 1/2, 1/2];
bq = [];                             Bq = [];


%fox = zeros(1, points, 'gpuArray');
fox = zeros(1, points);
for i = 1:points
    fox(i) = real(HFox(an, An, ap, Ap, bm, Bm ,bq, Bq, Xi(i)));
    %fox(i) = real(gpuArray(HFox(an, An, ap, Ap, bm, Bm ,bq, Bq, Xi(i))));
end

% compute CAP
CAP = preH .* fox;

end