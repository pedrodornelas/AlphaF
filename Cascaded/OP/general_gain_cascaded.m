function Gain = general_gain_cascaded(N, params, gammaBar)

% params is a matrix of params cascaded system
% params = [Channel 1: alpha, mu, ms, z, rc, Nc, Hl;
%           Channel 2: alpha, mu, ms, z, rc, Nc, Hl;
%           ...
%           Channel n: alpha, mu, ms, z, rc, Nc, Hl;]
% cascaded_gain = gain_channel1 * gain_channel2 * ... * gain_channeln

% N: number of channels
% alpha: non-linearity
% mu: number of multipaths
% ms: shadowing
% z: pointing error
% Ao: pointing error param
% rc: hat r (average power)
% Nc: number of points
% Hl: path loss
% gammaBar: SNR vector

channels = N;
columns = size(params, 2);

for c = 1:channels

    alpha = params(c,1);
    mu = params(c,2);
    ms = params(c,3);
    z = params(c,4);
    rc = params(c,5);
    Nc = params(c,6);
    Hl = params(c,7);

    Gain = ones(Nc, length(gammaBar));
    % Gain = 1;

    Ao = sqrt(gammaBar*(2+z^2))/(rc*z*Hl);
    for i = 1:length(gammaBar)
        % [c i]
        % random gains
        Hf = gainAF(alpha, mu, ms, rc, Nc, -1e-3); % Gain Alpha-F
        Hp = PointError(z, Ao(i), Nc);                % Gain Pointing Error
        % gain
        Gain(:,i) = Gain(:,i) .* (Hl(:) .* Hf(:) .* Hp(:));
        % line = size(Gain, 1);
        % col = size(Gain, 2);
    end

end

% Gain
% line = size(Gain,1)
% col = size(Gain,2)
Gain = Gain.';
% length(Gain,:)


end