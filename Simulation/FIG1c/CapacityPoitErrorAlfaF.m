clear all
clc

% Piece of code to calculate the Outage Probability
% Brito - 28/09/2022

set(0,'defaulttextinterpreter','latex');

% Número de amostras
Nc = 1e4;
bounds = [0 30];
points = 1e2;

% SNRs -- Amostragem dos valores observáveis
GBdB = linspace(0,30,15); % SNR em dB
gammaBar = 10.^(0.1*GBdB); % SNR linear

% Parâmetros da Distribuição Alfa F
alfa = [2, 5];
mu = 3;
ms = 5;
rc = 1;


% Perda de percurso
Hl = 1.00;

% Parametros da distribuição do erro de apontamento
z = [0.8, 1.5, 11];
%Ao = sqrt(gammaBar*(2+z(1)^2))/(rc*z(1)*Hl);

% % SNRs -- Amostragem dos valores observáveis
% GBdB = linspace(0,50,20); % SNR em dB
% gammaBar = 10.^(0.1*GBdB); % SNR linear
% rc = sqrt(gammaBar*(2+z^2))/(Ao*z*Hl);

% Inicialização dos vetores -- Prealocation
C = zeros(length(gammaBar));
Ao = zeros(length(gammaBar));

coloralfa = 'gb';

for Alfa = 1:length(alfa)
    for Z = 1:length(z)
        Ao = sqrt(gammaBar*(2+z(Z)^2))/(rc*z(Z)*Hl);
        for i = 1:length(gammaBar)
            i
            % Ganhos aleatórios
            Hf = gainAF(alfa(Alfa),mu,ms,rc,Nc,1e-3); % Alpha F
            Hp = PointError(z(Z),Ao(i),Nc); % Pointing error
            % Ganho total
            Gain = (Hl(:).*Hf(:).*Hp(:));

            % Amostral Ergodic Capacity 
            C(i) = mean(log2(1+Gain.^2)).';
        end
        [gammaBar_dB, P] = Capacity_asymptotic(alfa(Alfa), mu, ms, bounds, points, z(Z));
        [gammaBar_dB, Pb] = Capacity_analit(alfa(Alfa), mu, ms, bounds , points, z(Z));
        figure(1)
        plot(GBdB,C(:,1),'rx',...
                 gammaBar_dB,P,'k--',...
                 gammaBar_dB,Pb, coloralfa(Alfa),...
                 'linewidth', 1.2)

        hold on
    end
end

axis([min(GBdB) max(GBdB) 0 10])
legend('Simulated',"\alpha = 2",'','','','','','','','',"\alpha = 3.5", 'Asymptotic', 'Location', 'southwest')
%legend('Simulated','Asymptotic','\alpha=2','','','','','','','','', '\alpha = 5', 'Location', 'northwest')
%title('Capacity')
ylabel("Capacity")
xlabel("SNR")
grid on

%textbox com valores
dim = [0.15 0.25 0.2 0.2];
str = {"$\mu = 3$","$m_s = 5$"};
annotation('textbox',dim,'String',str,'FitBoxToText','on');
