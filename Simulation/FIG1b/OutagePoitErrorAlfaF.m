clear all
clc

% Piece of code to calculate the Outage Probability
% Brito - 28/09/2022


% SNRs -- Amostragem dos valores observáveis
GBdB = linspace(0,20,15); % SNR em dB
gammaBar = 10.^(0.1*GBdB); % SNR linear

% Limiares
GthdB = linspace(-5,10,2);
gth = 10.^(0.1*GthdB); % Limiar de SNR

% Parâmetros da Distribuição Alfa F
alfa = 3.5;
mu = 3;
ms = 5;
Nc = 1e4;
rc = 1;

% Perda de percurso
Hl = 1.00;

% Parametros da distribuição do erro de apontamento
z = 0.8;
Ao = sqrt(gammaBar*(2+z^2))/(rc*z*Hl);


% Inicialização dos vetores -- Prealocation
Pout = zeros(length(gammaBar),length(gth));

tic


for g = 1:length(gth)
    for i = 1:length(gammaBar)
        [g i]
        % Ganhos aleatórios
        Hf = gainAF(alfa,mu,ms,rc,Nc,1e-3); % Alpha F
        Hp = PointError(z,Ao(i),Nc); % Pointing error
        % Ganho total
        Gain = (Hl(:).*Hf(:).*Hp(:)).';
        
        [flagOP,~] = find(Gain.^2 <= gth(g));
        Pout(i,g) = sum(flagOP)/Nc;
        figure(1)
        subplot(2,1,1)
        title('OP x SNR (dB)')
        semilogy(GBdB,Pout(:,g),'-x'); hold on
        axis([min(GBdB) max(GBdB) min(Pout(:)) max(Pout(:))])
        subplot(2,1,2)
        title('OP x SNRth (dB)')
        semilogy(GthdB,Pout(i,:),'-x'); hold on
        axis([min(GthdB) max(GthdB) min(Pout(:)) max(Pout(:))])
    end
end

hold off
toc

%%
clc

% save CoisasQueHugerlesPediu

% L = 1 2 3
% gth = 5dB
% ms = 1.5 10 50 
% mu = 2.5 alpha = 2.5
[gammaBar_dB, Pb] = OP_analit(alfa, mu, ms, [1 gammaBar(end)], 100, z, gth(1));

clc

% % Parâmetros de entrada da MeijerG
Psi = mu/(ms-1);
m = 2; n = 2;
p = 3; q = 3;
a = [1-ms 1 1+z*z/alfa];
b = [mu z*z/alfa 0];
% PDF SNR
cdfG =@(gth,gb) z*z/alfa/gamma(mu)/gamma(ms).*...
          meijerG(a(1:n),a(n+1:end),...
                  b(1:m),b(m+1:end),...
                  Psi*(z*sqrt(gth./gb/(2+z^2))).^alfa);
figure(2)
semilogy(GBdB,Pout(:,1),'x',...
         GBdB,cdfG(gth(1),gammaBar),'-',...
         gammaBar_dB,Pb,'-.',...
         'linewidth',1.5)
axis([min(GBdB) max(GBdB) 1e-1 1e0])
grid on


%%
figure(2)
semilogy(GBdB,Cout(:),'-d',...
         GBdB,Cout(:),'-x',...
         GBdB,Cout(:),'-o')

