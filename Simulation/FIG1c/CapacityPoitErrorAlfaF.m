clear all
close all
clc

% Piece of code to calculate the Outage Probability
% Brito - 28/09/2022

set(0,'defaulttextinterpreter','latex');
p = [0,0,0; 0,0,0; 0,0,0];
q = [0,0,0; 0,0,0; 0,0,0];

% Número de amostras
Nc = 1e4;
bounds = [0 30];
points = 1e2;

% SNRs -- Amostragem dos valores observáveis
GBdB = linspace(0,30,15); % SNR em dB
gammaBar = 10.^(0.1*GBdB); % SNR linear

% Parâmetros da Distribuição Alfa F
alfa = 2.2;
mu = 2.1;
ms = [3.1, 10.5];
rc = 1;


% Perda de percurso
Hl = 1.00;

% Parametros da distribuição do erro de apontamento
z = [0.6, 1.1, 6.5];
%Ao = sqrt(gammaBar*(2+z(1)^2))/(rc*z(1)*Hl);

% % SNRs -- Amostragem dos valores observáveis
% GBdB = linspace(0,50,20); % SNR em dB
% gammaBar = 10.^(0.1*GBdB); % SNR linear
% rc = sqrt(gammaBar*(2+z^2))/(Ao*z*Hl);

% Inicialização dos vetores -- Prealocation
C = zeros(length(gammaBar));
Ao = zeros(length(gammaBar));

colorMS = 'gb';

for MS = 1:length(ms)
    for Z = 1:length(z)
        Ao = sqrt(gammaBar*(2+z(Z)^2))/(rc*z(Z)*Hl);
        for i = 1:length(gammaBar)
            Mzi = [MS Z i]
            % Ganhos aleatórios
            Hf = gainAF(alfa,mu,ms(MS),rc,Nc,1e-3); % Alpha F
            Hp = PointError(z(Z),Ao(i),Nc); % Pointing error
            % Ganho total
            Gain = (Hl(:).*Hf(:).*Hp(:));

            % Amostral Ergodic Capacity 
            C(i) = mean(log2(1+Gain.^2)).';
        end

        [gammaBar_dB, P] = Capacity_asymptotic(alfa, mu, ms(MS), bounds, points, z(Z));
        [gammaBar_dB, Pb] = Capacity_analit(alfa, mu, ms(MS), bounds , points, z(Z));
        figure(1)
 
        
        if MS == 1
        p(:, Z) = plot(GBdB,C(:,1),'rx',...
             gammaBar_dB,Pb, colorMS(MS),...
             gammaBar_dB,P,'k--',...
             'linewidth', 1.2 )
        end
        
        if MS == 2
        q(:, Z) = plot(GBdB,C(:,1),'rx',...
             gammaBar_dB,Pb, colorMS(MS),...
             gammaBar_dB,P,'k--',...
             'linewidth', 1.2 )
        end
        
        hold on
    end
end

axis([min(GBdB) max(GBdB) 0 10])
%legend('Simulated',"\alpha = 2",'','','','','','','','',"\alpha = 3.5", 'Asymptotic', 'Location', 'southwest')
legend([p(1) p(2) q(2) q(3)],'Simulated',"$m_s="+num2str(ms(1))+"$",'','','','','','','','', "$m_s="+num2str(ms(2))+"$","Asymptotic",'','Location','northwest')
set(legend, 'Interpreter', 'latex')
%title('Capacity')
ylabel("Capacity", 'FontSize', 14)
xlabel("SNR", 'FontSize', 14)
grid on

%textbox com valores
dim = [0.15 0.5 0.2 0.2];
str = {"$\mu="+num2str(mu)+"$","$\alpha="+num2str(alfa)+"$ (Fisher-Snedecor)"};
annotation('textbox',dim,'Interpreter','latex','String',str,'FitBoxToText','on');
