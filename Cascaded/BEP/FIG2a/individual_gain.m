function Gains = individual_gain(N, params, Nc, gammaBar)

% params is a matrix of params cascaded system
% params = [Channel 1: alpha, mu, ms, z, Ao, Hl;
%           Channel 2: alpha, mu, ms, z, Ao, Hl;
%           ...
%           Channel n: alpha, mu, ms, z, Ao, Hl;]
% cascaded_gain = gain_channel1 * gain_channel2 * ... * gain_channeln

% N: number of channels
% alpha: non-linearity
% mu: number of multipaths
% ms: shadowing
% z: pointing error
% Ao: pointing error param
% rc: hat r (average power)
% Nc: number of simulation points
% Hl: path loss
% gammaBar: SNR vector

channels = N;
columns = size(params, 2);

Gains = ones(Nc, length(gammaBar), N);

for c = 1:channels
    % c

    alpha = params(c,1);
    mu = params(c,2);
    ms = params(c,3);
    z = params(c,4);
    Ao = params(c,5);
    Hl = params(c,6);

    % Gain = 1;

    % Ao = sqrt(gammaBar(:, c) .* (2+z^2)) ./ (rc*z*Hl);
    rc = sqrt(1*(2+z^2))/(Ao*z*Hl);
    for i = 1:length(gammaBar)
        [c i]
        % random gains
        Hf = gainAF(alpha, mu, ms, rc, Nc, 1e-3); % Gain Alpha-F
        Hp = PointError(z, Ao, Nc);             % Gain Pointing Error
        % size(Gains(:,i,c))
        % size(Hl(:) .* Hf(:) .* Hp(:))
        Gains(:,i,c) = (Hl(:) .* Hf(:) .* Hp(:));
    end

end

end