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
alfa = [2.0, 3.5];
mu = 3;
ms = 1.3;
Nc = 1e4;
rc = 1;

% Perda de percurso
Hl = 1.00;

% Rrro de apontamento
z = [0.8, 1.5, 6.7];

% Inicialização dos vetores -- Prealocation
Pout = zeros(length(gammaBar));
Ao = zeros(length(gammaBar));

coloralfa = 'gb';

for j = 1:length(alfa)
    for k = 1:length(z)
        Ao = sqrt(gammaBar*(2+z(k)^2))/(rc*z(k)*Hl);
        for i = 1:length(gammaBar)
            i
            % Ganhos aleatórios
            Hf = gainAF(alfa(j),mu,ms,rc,Nc,1e-3); % Alpha F
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


axis([min(GBdB) max(GBdB) 1e-2 1])
legend('Simulated',"$\alpha = 2$",'','','','','','','','',"$\alpha = 3.5$", 'Asymptotic', 'Location', 'southwest')
%title('Outage Probability')
set(legend, 'Interpreter', 'latex')
ylabel("OP", 'FontSize', 14)
xlabel("SNR", 'FontSize', 14)
grid on

%textbox com valores
dim = [0.15 0.25 0.2 0.2];
str = {"$\mu = 3$","$m_s = 1.3$"};
annotation('textbox',dim,,'interpreter','latex','String',str,'FitBoxToText','on');

%GBdB,cdfG(gth(1),gammaBar),'-',