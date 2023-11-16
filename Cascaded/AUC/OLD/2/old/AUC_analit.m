function [gammaBar_dB, auc] = AUC_analit(N, alpha, mu, ms, bounds, points, z, u)
% TODO: doc func
% points - number of points

%lambda = 1;

% generate independent variable vector
L = bounds(1);
U = bounds(2);
gammaBar_dB = linspace(L, U, points);
%gammaBar = gpuArray.linspace(L, U, points);

gammaBar = db2pow(gammaBar_dB);

onesN = ones(1, N);
auc = ones(1, length(gammaBar));

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


        % m = 2N; n = N+1; p = 2N+1; q = 2N;
        an = [(1-l), (1-ms)*onesN];              An = [1/2, (1/alpha)*onesN];   
        ap = [((z.^2)/alpha) + 1];               Ap = [(1/alpha)*onesN];
        bm = [mu, ((z.^2)/alpha)];               Bm = [(1/alpha)*onesN, (1/alpha)*onesN];
        bq = [];                                 Bq = [];
    
    
        %fox = zeros(1, points, 'gpuArray');
        fox = zeros(1, points);
        for i = 1:points
            fox(i) = real(HFox(an, An, ap, Ap, bm, Bm ,bq, Bq, Xi(i)));
            %fox(i) = real(gpuArray(HFox(an, An, ap, Ap, bm, Bm ,bq, Bq, Xi(i))));
        end

        % compute AUC
        auc = auc - preH .* fox;
    end
end

% DEBUG
% semilogy(gammaBar_dB, Pb)
end