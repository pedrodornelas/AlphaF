function [gammaBar_dB, auc_asymp] = AUC_asymptotic(N, alpha, mu, ms, bounds, points, z, u)
% Function to implement the asymptotic expression 
% for the alpha-F fading and AWGGN distributions.
%
% INPUTS:
% alpha - positive power parameter
% beta - 
% mu - number of multipath clusters
% ms - shadowing parameter
% points - number of points
%
% OUTPUTS:
% P - points that form the asymptote

% generate independent variable vector
L = bounds(1);
U = bounds(2);
gammaBar_dB = linspace(L, U, points);
%gammaBar = gpuArray.linspace(L, U, points);

gammaBar = db2pow(gammaBar_dB);
auc_asymp = ones(1, length(gammaBar));

for k = 0:(u-1)
    for l = 0:k
        Xi = ones(1, length(gammaBar));
        preH = 1;
        for j = 1:N
            Psi = (mu(j)/(ms-1)) ^ (1/alpha);
            Xi = Xi .* (Psi .* (z(j) ./ (sqrt(gammaBar .* (z(j)^2+2)))));

            % precomputations
            preH = preH * (z(j)^2/(alpha*gamma(mu(j))*gamma(ms)));
        end

        Xi = Xi .* sqrt(2);
        preH = preH * (1/((2^(k+u+1)) * factorial(l))) * nchoosek((k+u-1), (k-l));


        onesN = ones(1, N);

        % m = 2N; n = N+1; p = 2N+1; q = 2N;
        an = [(1-l), (1-ms)*onesN];              An = [1/2, (1/alpha)*onesN];   
        ap = [((z.^2)/alpha) + 1];               Ap = [(1/alpha)*onesN];
        bm = [mu, ((z.^2)/alpha)];               Bm = [(1/alpha)*onesN, (1/alpha)*onesN];
        bq = [];                                 Bq = [];


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

        P = (preH * phiU * (Xi .^ U)) ./ Bc;
        auc_asymp = auc_asymp - P;

    end
end

end