function [AUC] = AUC_asymptotic(N, params, u, gammaBar)
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

AUC = 1;

for k = 0:(u-1)
    for l = 0:k

        Xi = 1;
        preH = 1;
        for c = 1:N
            alpha = params(c,1);
            mu = params(c,2);
            ms = params(c,3);
            z = params(c,4);

            Psi = (mu/(ms-1)) ^ (1/alpha);
            Xi = Xi .* ((Psi .* z) ./ sqrt( gammaBar(:, c) .* (z^2+2)));

            % precomputations
            preH = preH * (z^2/(alpha * gamma(mu) * gamma(ms)));
        end

        Xi = (Xi .* sqrt(2))';
        preH = preH * (1/((2^(k+u+1)) * factorial(l))) * nchoosek((k+u-1), (k-l));

        alphas = params(1:N, 1).';
        mus = params(1:N, 2).';
        mss = params(1:N, 3).';
        zs = params(1:N, 4).';

        % m = 2N; n = N+1; p = 2N+1; q = 2N;
        an = [(1-l), (1-mss)];              An = [1/2, (1./alphas)];   
        ap = [((zs.^2)./alphas) + 1];       Ap = [(1./alphas)];
        bm = [mus, ((zs.^2)./alphas)];      Bm = [(1./alphas), (1./alphas)];
        bq = [];                            Bq = [];


        div = bm./Bm;
        [U, idx] = min(div);

        mins = find(div == U);

        Bc = 1;
        Bc = Bm(idx);

        % gamma 1
        gamma1 = 1;
        for j = 1:length(bm)
            % j not is member of mins
            if ismember(j,mins) == 0
                gamma1 = gamma1*gamma(bm(j) - U*Bm(j));
            else
                % Bc = Bc*Bm(j);
            end
        end

        gamma2 = 1;
        for j = 1:length(an)
            gamma2 = gamma2*gamma(1 - an(j) + U*An(j));
        end

        gamma3 = 1;
        for j = 1:length(bq)
            gamma3 = gamma3*gamma(1 - bq(j) + U*Bq(j));
        end

        gamma4 = 1;
        for j = 1:length(ap)
            gamma4 = gamma4*gamma(ap(j) - U*Ap(j));
        end

        phiU = (gamma1*gamma2) / (gamma3*gamma4);

        P = (preH .* phiU .* (Xi .^ U)) ./ Bc;
        AUC = AUC - P;

    end
end

end