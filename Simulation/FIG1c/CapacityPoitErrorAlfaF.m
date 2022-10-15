clear all
clc

% Piece of code to calculate the Outage Probability
% Brito - 28/09/2022

% Número de amostras
Nc = 1e4;
bounds = [0 30];
points = 1e5;

% SNRs -- Amostragem dos valores observáveis
GBdB = linspace(0,30,15); % SNR em dB
gammaBar = 10.^(0.1*GBdB); % SNR linear

% Parâmetros da Distribuição Alfa F
alfa = 3.5;
mu = 3;
ms = 5;
rc = 1;

% Perda de percurso
Hl = 1.00;

% Parametros da distribuição do erro de apontamento
z = [0.8, 1.5, 6.7];
%Ao = sqrt(gammaBar*(2+z(1)^2))/(rc*z(1)*Hl);

% % SNRs -- Amostragem dos valores observáveis
% GBdB = linspace(0,50,20); % SNR em dB
% gammaBar = 10.^(0.1*GBdB); % SNR linear
% rc = sqrt(gammaBar*(2+z^2))/(Ao*z*Hl);

% Inicialização dos vetores -- Prealocation
C = zeros(length(gammaBar),1);

tic
for Z = 1:length(z)
    Ao = sqrt(gammaBar*(2+z(Z)^2))/(rc*z(Z)*Hl);
    for i = 1:length(gammaBar)
        i
        % Ganhos aleatórios
        Hf = gainAF(alfa,mu,ms,rc,Nc,1e-3); % Alpha F
        Hp = PointError(z(Z),Ao(i),Nc); % Pointing error
        % Ganho total
        Gain = (Hl(:).*Hf(:).*Hp(:));

        % Amostral Ergodic Capacity 
        C(i) = mean(log2(1+Gain.^2)).';
        figure(1)
        title('C x SNR (dB)')
        semilogy(GBdB,C,'-x'); hold on
        axis([min(GBdB) max(GBdB) min(C(:)) max(C(:))])
        
    end
    [gammaBar_dB, P] = Capacity_asymptotic(alfa, mu, ms, bounds, points, z(Z));
    [gammaBar_dB, Pb] = Capacity_analit(alfa, mu, ms, bounds , points, z(Z));
    figure(2)
    plot(GBdB,C(:),'o',...
             gammaBar_dB,P,'r',...
             gammaBar_dB,Pb,'b')
    hold on
end

hold off
toc

%%

% [gammaBar_dB, P] = Capacity_asymptotic(alfa, mu, ms, bounds, points, z(1));
% [gammaBar_dB, Pb] = Capacity_analit(alfa, mu, ms, bounds , points, z(1));
% figure(2)
% plot(GBdB,C(:),'o',...
%          gammaBar_dB,P,'r',...
%          gammaBar_dB,Pb,'b')
