clear all
close all
clc

% Piece of code to calculate the Outage Probability
% Pedro Henrique Dornelas Almeida - 05/05/2023

set(0,'defaulttextinterpreter','latex');
set(groot,'defaultAxesTickLabelInterpreter','latex');

% SNRs -- Amostragem dos valores observáveis
L = 0;    %db
U = 30;   %db
points = 1e2;
bounds = [L U]; %dB gammaBar limits

GBdB = linspace(L, U, 15); % SNR em dB
gammaBar = 10.^(0.1*GBdB); % SNR linear

% Parâmetros da Distribuição Alpha F - Cascaded
N = 1;              % nº estações relay
alpha = [2, 2.5];
mu = 0.3;
ms = 90;

% Parâmetro AUC
u = 2;

% Erro de apontamento
% z = [0.6, 1.0, 1.5, 13];
% z = [0.6, 0.8; 1.0, 1.1; 3, 3.3; 14, 15];
% z = [0.6, 0.8; 1.0, 1.1; 1.5, 1.6; 3, 3.3; 14, 15];
z = [0.6, 0.8; 1.0, 1.1; 1.5, 1.6; 14, 15];

for i=1:length(alpha)
    if ms <= (4/alpha(i))
        fprintf('ms > 4/alpha not met\n')
        error
    end
end

% Parâmetros da simulação
Nc = 1e4; % Número de pontos
rc = 1;   % 

% Perda de percurso
Hl = 1.00;

% Inicialização dos vetores -- Prealocation
Pout = zeros(length(gammaBar));
Ao = zeros(length(gammaBar));

colorz = 'brgmp';

for w = 1:length(z)
        
    [gammaBar_dB, auc] = AUC_analit(N, alpha(N), mu, ms, bounds, points, z(w, N), u);
    [gammaBar_dB, auc_asymp] = AUC_asymptotic(N, alpha(N), mu, ms, bounds, points, z(w, N), u);

    figure(1)
    % plot(GBdB, Pout(:,1),'rx',...
    %         gammaBar_dB, auc, colorz(j),...
    %         gammaBar_dB, P,'w--',...
    %         'linewidth',1.2)
    % str = colorz(j)+'o'
        
    % if w == 4
    %     color = 'm--'
    % else
    %     color = colorz(w)
    % end
    plot(gammaBar_dB, auc, colorz(w),...
         gammaBar_dB, auc_asymp,'k--',...
         'linewidth', 1.2)

    % plot(gammaBar_dB, auc, colorz(w),...
    %      'linewidth', 1.2)    
    
    hold on
end

axis([min(GBdB) max(GBdB) 0.5 1])
tam_fonte = 11;
legend('FontSize', tam_fonte)
% legend("$z_1="+num2str(z(1,1))+", z_2="+num2str(z(1,2))+"$", '', "$z_1="+num2str(z(2,1))+", z_2="+num2str(z(2,2))+"$", '', "$z_1="+num2str(z(3,1))+", z_2="+num2str(z(3,2))+"$", '', "$z_1="+num2str(z(4,1))+", z_2="+num2str(z(4,2))+"$", 'Asymptotic', 'Location', 'southwest')
legend("$z_1="+num2str(z(1,1))+", z_2="+num2str(z(1,2))+"$", '',"$z_1="+num2str(z(2,1))+", z_2="+num2str(z(2,2))+"$", '',"$z_1="+num2str(z(3,1))+", z_2="+num2str(z(3,2))+"$", '', "Nakagami-m", 'Asymptotic', 'Location', 'southeast')
ax = gca;
ax.FontSize = 13;

%title('Outage Probability')
set(legend, 'Interpreter', 'latex')
ylabel("AUC", 'FontSize', 14)
yticks([0.5:0.05:1])
xlabel("SNR (dB)", 'FontSize', 14)
grid on

%textbox com valores
dim = [0.75 0.40 0.2 0.2];
str = {"$\alpha_1 ="+num2str(alpha(1))+"$", "$\alpha_2 ="+num2str(alpha(2))+"$","$\mu ="+num2str(mu)+"$", "$m ="+num2str(ms)+"$", "$u ="+num2str(u)+"$"};
annotation('textbox',dim,'interpreter','latex','String',str,'FitBoxToText','on', 'FontSize', tam_fonte);