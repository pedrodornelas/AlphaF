function [sep,bep] = channel(M, N, params, Nc, gammaBar)
% params is a matrix of params cascaded system
% params = [Channel 1: alpha, mu, ms, z, Ao, Hl;
%           Channel 2: alpha, mu, ms, z, Ao, Hl;
%           ...
%           Channel n: alpha, mu, ms, z, Ao, Hl;]
% cascaded_gain = gain_channel1 * gain_channel2 * ... * gain_channeln

% M: constelation order
% N: number of channels
% alpha: non-linearity
% mu: number of multipaths
% ms: shadowing
% z: pointing error
% Ao: pointing error param
% rc: hat r (average power)
% Nc: number of simulation points (bits)
% Hl: path loss
% gammaBar: SNR matriz per channel
% gammaBar = [          , channel1, channel2, ..., channeln;
%              gammaBar1,     y1.1,     y1.2, ...,     y1.n;
%                    ...,      ...,      ..., ...,      ...;
%             gammaBar15,    y15.1,    y15.2, ...,    y15.n;]

k = log2(M); % Numero de bits por simbolo
if mod(Nc,log2(M)) ~= 0 % Correcao no numero de bits
    Nc = Nc+log2(M)-mod(Nc,log2(M));
end

coded = randi([0 1],Nc,1);
dados = coded; % Geracao da sequencia binaria
dados_s2p = reshape(dados,Nc/k,k); % Conversao Serie-Paralelo
dados_dec = bi2de(dados_s2p); % Conversao binario-decimal para futura
                              % correlacao coms os simbolos da
                              % constelacao

s = qammod(0:1:M-1,M,'gray','UnitAveragePower',true); % Geracao dos simbolos da constelacaoo
                                                      % utilizando mapeamento Gray

qtd_dados = length(dados_dec);


% ga = 10.^(SNR_dB/10);
n = zeros(length(gammaBar), qtd_dados);
for i = 1:length(gammaBar)
    % gammaBar(i, 1)
    n(i, :) = 1j*normrnd(0,1/sqrt(2*gammaBar(i, N)),[1 qtd_dados])+...
                 normrnd(0,1/sqrt(2*gammaBar(i, N)),[1 qtd_dados]);
    
end

gain_channels = individual_gain(N, params, qtd_dados, gammaBar);
gain_total = cascaded_gain(gain_channels);


% rc = sqrt(1*(2+z^2))/(Ao*z*Hl);
% Hf = gainAF(alfa,mu,ms,rc,Nc,-1e-3);
% Hp = PointError(z,Ao,Nc);
% Gain = (Hl(:).*Hf(:).*Hp(:)).';

Gain = gain_total(:, :, N).';

r = Gain.*s(dados_dec+1) + n;

bep = zeros(1, length(gammaBar));
for i = 1:length(gammaBar)
    r_i = r(i,:);
    demod = qamdemod(r_i,M,'gray');

    dados_demod = de2bi(demod,k);
    dados_p2s = dados_demod(:);

    [~,bep(i)] = biterr(coded,dados_p2s);
end

sep = bep.*log2(M);

end