function [OP] = OP_analit(L, N, params, gamma_th, gammaBar)
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
    
    Xi = (Xi .* sqrt(gamma_th)); % ./ sqrt(gammaBar);

    preH = (preH ^ L) / 2;

    % Python function for multivariate Fox H
    py.importlib.reload(py.importlib.import_module('multiFoxH'));
    pyModule = py.importlib.import_module('multiFoxH');
    H = double(pyModule.parseArgsFromMatlab(params, N, L, Xi));
    
    % compute OP
    OP = preH .* H;
end