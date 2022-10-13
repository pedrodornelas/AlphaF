function [ser,ber] = channel(M,SNR_dB,N,alfa,mu,ms,z,Ao,Hl)

ga = 10.^(SNR_dB/10);

k = log2(M); % N�mero de bits por s�mbolo
if mod(N,log2(M)) ~= 0 % Corre��o no n�mero de bits
    N = N+log2(M)-mod(N,log2(M));
end

coded = randi([0 1],N,1);
dados = coded; % Gera��o da sequ�ncia bin�ria
dados_s2p = reshape(dados,N/k,k); % Convers�o S�rie-Paralelo
dados_dec = bi2de(dados_s2p); % Convers�o binario-decimal para futura
                              % correla��o coms os s�mbolos da
                              % constela��o

s = qammod(0:1:M-1,M,'gray','UnitAveragePower',true); % Gera��o dos s�mbolos da constela��o
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