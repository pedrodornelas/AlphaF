
clc
clear all
close all

% Piece of code to calculate the AUC
% Pedro Henrique Dornelas Almeida - 10/11/2023

functions_path = "functions";
addpath(functions_path);

set(0,'defaulttextinterpreter','latex');
set(groot,'defaultAxesTickLabelInterpreter','latex');

% SNRs -- Amostragem dos valores teoricos
L = 0;    %db
U = 50;   %db
points = 1e2;
bounds = [L U]; %dB gammaBar limits

analit_gammaBar_dB = linspace(L, U, points);
analit_gammaBar = 10.^(0.1.*analit_gammaBar_dB);
simu_gammaBar_dB = linspace(L, U, 15); % SNR em dB
simu_gammaBar = 10.^(0.1.*simu_gammaBar_dB); % SNR linear

% Parâmetros da Distribuição Alpha-F - Cascaded
% Parâmetros para obter Rayleigh
N = 1;    % nº estações relay
if N == 1
    ASY = 1;
else
    ASY = 0;
end

% curva 1: nakagami-m
% curva 2: fisher
% curva 3: weibull
curvas = 3;
alpha = [2, 2, 3.1];
mu = [0.6, 0.61; 
      0.8, 0.81; 
      1.0, 1.01;];
ms = [3.5];                   % shadowing -> infty
z = [0.7, 0.72;
     8.0, 9.00;];            % pointing error -> infty
% z = [1.2, 1.21];              % pointing error -> infty

% Parâmetro AUC
u = 4;
N_auc = 50;
N_roc = 5e2;

% gammaBar por canal
analit_gammaBar_c = db2pow(1) * (ones( length(analit_gammaBar) , max(N)));
analit_gammaBar_c(:, 1) = analit_gammaBar; % variar só o primeiro canal, os demais mantemos em 1dB
simu_gammaBar_c = db2pow(1) * (ones( length(simu_gammaBar) , max(N)));
simu_gammaBar_c(:, 1) = simu_gammaBar;     % variar só o primeiro canal, os demais mantemos em 1dB


for i=1:length(alpha)
    if ms <= (4/alpha(i))
        fprintf('ms > 4/alpha not met\n')
        error
    end
end

% Parametros da simulacao
Nc = 1e4;  % numero de pontos para simulacao de monte carlo
Ao = 0.8;
rc = 1;    % hat r
Hl = 1.00; % loss path

colorz = 'brgmp';
tic;
h = [];
cont = 1;

for k = 1:size(z, 1)
    for c = 1:curvas
        c

        simulation_params = [alpha(c), mu(c, 1), ms, z(k, 1), rc, Hl;
                            alpha(c), mu(c, 2), ms, z(k, 2), rc, Hl;];

        analit_params = [alpha(c), mu(c, 1), ms, z(k, 1);
                        alpha(c), mu(c, 2), ms, z(k, 2);];

        sim_AUC = AUC_simulation(Nc, N, simulation_params, simu_gammaBar_c, u, N_auc, N_roc);
        [auc] = AUC_analit(N, analit_params, u, analit_gammaBar_c);  % teorica
        [asy_auc] = AUC_asymptotic(N, analit_params, u, analit_gammaBar_c);

        figure(1)

        h(cont) = plot(analit_gammaBar_dB, auc, colorz(c),'linewidth', 1.2); hold on;
        cont = cont+1;
        if ASY == 1
            h(cont) = plot(analit_gammaBar_dB, asy_auc, ['k--'],'linewidth', 0.8); hold on;
            cont = cont+1;
        end
        h(cont) = plot(simu_gammaBar_dB, sim_AUC, 'rx','linewidth', 1.2); hold on;
        cont = cont+1;         
        
        hold on
    end
end

execution_time = toc;
format_time = executionTime(execution_time);
disp(['Execution time: ' format_time]);

axis([min(simu_gammaBar_dB) max(simu_gammaBar_dB) 0.5 1])
tam_fonte = 11;
legend('FontSize', tam_fonte)
curvas_str = ["Fisher-Snedecor ${\cal F}$ ($\mu="+num2str(mu(1,1))+"$)", "Fisher-Snedecor ${\cal F}$ ($\mu="+num2str(mu(2,1))+"$)", "Weibull-${\cal F}$ ($\alpha="+num2str(alpha(3))+"$)"];
% legend(curvas_str(1), '','','', curvas_str(2), '','','',curvas_str(3), 'Asymptotic', 'Location', 'southeast')
% legend(curvas_str(1), '', '', curvas_str(2), '', '',curvas_str(3), 'Asymptotic', 'Simulated', 'Location', 'southeast')
legend("Location", "southeast");
if ASY == 1
    legend([h(1), h(4), h(7), h(8), h(9)], {curvas_str(1), curvas_str(2), curvas_str(3), "Asymptotic", "Simulated"});
else
    legend([h(1), h(3), h(5), h(6)], {curvas_str(1), curvas_str(2), curvas_str(3), "Simulated"});
end
ax = gca;
ax.FontSize = 13;

%title('Outage Probability')
set(legend, 'Interpreter', 'latex')
ylabel("Average AUC", 'FontSize', 14)
yticks([0.5:0.05:1])
xlabel("SNR (dB)", 'FontSize', 14)
grid on

%textbox com valores
dim = [0.755 0.37 0.2 0.2];
str = {"$z ="+num2str(z(1,1))+"$", "$N ="+num2str(N(1))+"$", "$m ="+num2str(ms(1))+"$", "$u ="+num2str(u)+"$"};
annotation('textbox',dim,'interpreter','latex','String',str,'FitBoxToText','on', 'FontSize', tam_fonte, "BackgroundColor", "white");