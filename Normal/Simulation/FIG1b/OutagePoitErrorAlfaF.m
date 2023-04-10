clear all
close all
clc

% Piece of code to calculate the Outage Probability
% Brito - 28/09/202

set(0,'defaulttextinterpreter','latex');

% SNRs -- Amostragem dos valores observáveis
L = 0;    %db
U = 30;   %db
points = 1e2;
bounds = [L U]; %dB gammaBar limits

GBdB = linspace(L, U, 15); % SNR em dB
gammaBar = 10.^(0.1*GBdB); % SNR linear

% Limiares
% GthdB = linspace(-5,10,2);
% gth = 10.^(0.1*GthdB); % Limiar de SNR
gamma_th_dB = 5; %in dB
gamma_th = db2pow(gamma_th_dB);

% Parâmetros da Distribuição Alfa F
alfa = [2.0, 5.7];
mu = 2.7;
ms = 1.3;
Nc = 1e4;
rc = 1;

% Perda de percurso
Hl = 1.00;

% Rrro de apontamento
z = [0.6, 1.1, 6.5];

% Inicialização dos vetores -- Prealocation
Pout = zeros(length(gammaBar));
Ao = zeros(length(gammaBar));

coloralfa = 'gb';

for j = 1:length(alfa)
    for k = 1:length(z)
        Ao = sqrt(gammaBar*(2+z(k)^2))/(rc*z(k)*Hl);
        for i = 1:length(gammaBar)
            [j k i]
            % Ganhos aleatórios
            Hf = gainAF(alfa(j),mu,ms,rc,Nc,-1e-3); % Alpha F
            Hp = PointError(z(k),Ao(i),Nc); % Pointing error
            % Ganho total
            Gain = (Hl(:).*Hf(:).*Hp(:)).';
            
            [flagOP,~] = find(Gain.^2 <= gamma_th);
            Pout(i) = sum(flagOP)/Nc;
        end

        %%

        [gammaBar_dB, Pb] = OP_analit(alfa(j), mu, ms, bounds, points, z(k), gamma_th);
        [gammaBar_dB, P] = OP_asymptotic(alfa(j), mu, ms, bounds, points, z(k), gamma_th);

        figure(1)
        semilogy(GBdB, Pout(:,1),'rx',...
                gammaBar_dB, Pb, coloralfa(j),...
                gammaBar_dB, P,'k--',...
                'linewidth',1.2)
        hold on
    end
end


axis([min(GBdB) max(GBdB) 1e-5 1])
legend('Simulated',"$\alpha="+num2str(alfa(1))+"$ (Fisher-Snedecor)",'','','','','','','','',"$\alpha="+num2str(alfa(2))+"$", 'Asymptotic', 'Location', 'southwest')
%title('Outage Probability')
set(legend, 'Interpreter', 'latex')
ylabel("OP", 'FontSize', 14)
yticks([1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1e0])
xlabel("SNR", 'FontSize', 14)
grid on

%textbox com valores
dim = [0.15 0.25 0.2 0.2];
str = {"$\mu ="+num2str(mu)+"$" , "$m_s ="+num2str(ms)+"$", "$\gamma_{th}="+num2str(gamma_th_dB)+"$ dB"};
annotation('textbox',dim,'interpreter','latex','String',str,'FitBoxToText','on');

%GBdB,cdfG(gth(1),gammaBar),'-',