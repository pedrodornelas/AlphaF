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
L = -10;    %db
U = 30;   %db
points = 1e2;
bounds = [L U]; %dB gammaBar limits

analit_gammaBar_dB = linspace(L, U, points);
analit_gammaBar = 10.^(0.1.*analit_gammaBar_dB);
simu_gammaBar_dB = linspace(L, U, 15); % SNR em dB
simu_gammaBar = 10.^(0.1.*simu_gammaBar_dB); % SNR linear

% Parâmetros da Distribuição Alpha-F - Cascaded
% Parâmetros para obter Rayleigh
N = 1;              % nº estações relay
alpha = 2;
mu = 1;
ms = 80;            % shadowing -> infty
z = [8,15];         % pointing error -> infty

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
Nc = 5e3;  % numero de pontos para simulacao de monte carlo
rc = 1;    % hat r
Hl = 1.00; % loss path

simulation_params = [alpha, mu, ms, z(1), rc, Hl;];

sim_AUC = AUC_simulation_new(Nc, N, simulation_params, simu_gammaBar_c, u, N_auc, N_roc);

colorz = 'brgmp';

for w = 1:length(z)
    analit_params = [alpha, mu, ms, z(w);];
    
    [auc] = AUC_analit(N, analit_params, u, analit_gammaBar_c); % teorica

    figure(1)

    plot(analit_gammaBar_dB, auc, colorz(w),...
         simu_gammaBar_dB, sim_AUC, 'rx',...  
         'linewidth', 1.2)    
    
    hold on
end

axis([min(simu_gammaBar_dB) max(simu_gammaBar_dB) 0.5 1])
tam_fonte = 11;
legend('FontSize', tam_fonte);
legend("Rayleigh", "Simulated", "Location", "southeast");
ax = gca;
ax.FontSize = 13;

%title('Outage Probability')
set(legend, 'Interpreter', 'latex')
ylabel("Average AUC", 'FontSize', 14)
yticks([0.5:0.05:1])
xlabel("SNR (dB)", 'FontSize', 14)
grid on