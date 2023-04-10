function [ser,ber] = channel(M,SNR_dB,N,alfa,mu,ms,z,Ao,Hl)

ga = 10.^(SNR_dB/10);

k = log2(M); % Número de bits por símbolo
if mod(N,log2(M)) ~= 0 % Correção no número de bits
    N = N+log2(M)-mod(N,log2(M));
end

coded = randi([0 1],N,1);
dados = coded; % Geração da sequência binária
dados_s2p = reshape(dados,N/k,k); % Conversão Série-Paralelo
dados_dec = bi2de(dados_s2p); % Conversão binario-decimal para futura
                              % correlação coms os símbolos da
                              % constelação

s = qammod(0:1:M-1,M,'gray','UnitAveragePower',true); % Geração dos símbolos da constelação
                                                      % utilizando mapeamento Gray

Nc = length(dados_dec);

n = 1j*normrnd(0,1/sqrt(2*ga),[1 Nc])+...
       normrnd(0,1/sqrt(2*ga),[1 Nc]);

rc = sqrt(1*(2+z^2))/(Ao*z*Hl);
Hf = gainAF(alfa,mu,ms,rc,Nc,1e-3);
Hp = PointError(z,Ao,Nc);


Gain = (Hl(:).*Hf(:).*Hp(:)).';

r = s(dados_dec+1) + n./Gain;

demod = qamdemod(r,M,'gray');

dados_demod = de2bi(demod,k);
dados_p2s = dados_demod(:);

[~,ber] = biterr(coded,dados_p2s);
ser = length(find(demod(:)~=dados_dec(:)))/N;

end