clear all
close all
clc

% Piece of code to calculate the Outage Probability
% Pedro Henrique Dornelas Almeida - 05/05/2023

set(0,'defaulttextinterpreter','latex');
set(groot,'defaultAxesTickLabelInterpreter','latex');

% SNRs -- Amostragem dos valores observáveis
L = 0;    %db
U = 70;   %db
% points = 1e2;
points = 100;
bounds = [L U]; %dB gammaBar limits

analit_gammaBar_dB = linspace(L, U, points);
analit_gammaBar = 10.^(0.1.*analit_gammaBar_dB);
simu_gammaBar_dB = linspace(L, U, 15); % SNR em dB
simu_gammaBar = 10.^(0.1.*simu_gammaBar_dB); % SNR linear

% Parâmetros da Distribuição Alpha F - Cascaded
N = [2];              % nº estações relay

% curva 1: nakagami-m
% curva 2: fisher
% curva 3: weibull
curvas = 3;
alpha = [2, 2, 3.3];
mu = [0.3, 0.31 ; 0.6, 0.61 ; 1, 1.05];
ms = [3.5, 3.5, 3.5];

% Parâmetro AUC
u = 2;

analit_gammaBar_c = ones( length(analit_gammaBar) , max(N));
analit_gammaBar_c(:, 1) = analit_gammaBar; % variar só do primeiro canal...
simu_gammaBar_c = ones( length(simu_gammaBar) , max(N));
simu_gammaBar_c(:, 1) = simu_gammaBar; % variar só do primeiro canal...

% Erro de apontamento
% z = [0.6, 1.0, 1.5, 13];
% z = [0.6, 0.8; 1.0, 1.1; 3, 3.3; 14, 15];
% z = [0.6, 0.8; 1.0, 1.1; 1.5, 1.6; 3, 3.3; 14, 15];
z = [0.7, 0.77 ; 0.7, 0.77];

for i=1:length(alpha)
    if ms(i) <= (4/alpha(i))
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

for c = 1:curvas
    for w = 1:length(N)
        
        n = N(w);
        [gammaBar_dB, auc]       =     AUC_analit(N(w), alpha(c), mu(c, 1:n), ms(c), bounds, points, z(w, 1:n), u);
        [gammaBar_dB, auc_asymp] = AUC_asymptotic(N(w), alpha(c), mu(c, 1:n), ms(c), bounds, points, z(w, 1:n), u);

        figure(1)
        
        plot(gammaBar_dB, auc, colorz(c),...
             gammaBar_dB, auc_asymp,'k--',...
             'linewidth', 1.2)  
        
        hold on
    end
end

axis([min(GBdB) max(GBdB) 0.5 1])
tam_fonte = 11;
legend('FontSize', tam_fonte)
curvas_str = ["Nakagami-${\cal F}$ ($\mu="+num2str(mu(1,1))+"$)", "Fisher-Snedecor ${\cal F}$ ($\mu="+num2str(mu(2,1))+"$)", "Weibull-${\cal F}$ ($\alpha="+num2str(alpha(3))+"$)"];
% legend(curvas_str(1), '','','', curvas_str(2), '','','',curvas_str(3), 'Asymptotic', 'Location', 'southeast')
legend(curvas_str(1), '', curvas_str(2), '',curvas_str(3), 'Simulated', 'Location', 'southeast')
ax = gca;
ax.FontSize = 13;

%title('Outage Probability')
set(legend, 'Interpreter', 'latex')
ylabel("Average AUC", 'FontSize', 14)
yticks([0.5:0.05:1])
xlabel("SNR (dB)", 'FontSize', 14)
grid on

%textbox com valores
dim = [0.755 0.35 0.2 0.2];
str = {"$z ="+num2str(z(1,1))+"$", "$N ="+num2str(N(1))+"$", "$m ="+num2str(ms(1))+"$", "$u ="+num2str(u)+"$"};
annotation('textbox',dim,'interpreter','latex','String',str,'FitBoxToText','on', 'FontSize', tam_fonte, "BackgroundColor", "white");