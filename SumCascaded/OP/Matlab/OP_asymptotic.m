function [OP] = OP_asymptotic(L, N, params, gamma_th, gammaBar)
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

    sum_channels = L;
    channels = N;
    points = length(gammaBar);

    Xi = 1;
    preH = 1;
    for c = 1:N
        % c

        alpha = params(c,1);
        mu = params(c,2);
        ms = params(c,3);
        z = params(c,4);
    
        Psi = (mu/(ms-1)) ^ (1/alpha);
        Xi = Xi .* ((Psi .* z) ./ sqrt( gammaBar(:, c) .* (z^2+2)));
    
        % precomputations
        preH = preH * (z^2/(alpha * gamma(mu) * gamma(ms)));
    end

    Xi = (Xi .* sqrt(gamma_th))'; % ./ sqrt(gammaBar);

    preH = (preH ^ L) / 2;

    alphas = params(1:N, 1).';
    mus = params(1:N, 2).';
    mss = params(1:N, 3).';
    zs = params(1:N, 4).';

    an = [1, (1 - mss)];                 An = [1, (1./alphas)];   
    ap = [(((zs.^2)./alphas) + 1)];      Ap = [(1./alphas)];
    bm = [mus, ((zs.^2)./alphas)];       Bm = [(1./alphas), (1./alphas)];
    bq = [];                             Bq = [];

    div = bm./Bm;
    [U, idx] = min(div);

    mins = find(div == U);

    Bc = 1;
    Bc = Bm(idx) ^ L;

    % gamma 1
    gamma1 = 1;
    for j = 1:length(bm)
        % j is not member of mins
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

    phiU = ((gamma1*gamma2) / (gamma3*gamma4)) ^ L;

    gamma_sum1 = gamma((U/2) * L);
    gamma_sum2 = gamma((U*L));
    gamma_sum3 = gamma(1 + ((U/2) * L));
    psiU = gamma_sum1 / (gamma_sum2*gamma_sum3);


    OP = (preH * psiU * phiU * (Xi .^ (U*L))) ./ Bc;

end