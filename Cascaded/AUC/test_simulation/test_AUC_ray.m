clear all
close all
clc

% Piece of code to calculate the AUC
% Pedro Henrique Dornelas Almeida - 05/05/2023

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

% Parâmetros da Distribuição Alpha F - Cascaded
N = 1;              % nº estações relay
alpha = 2;
mu = 1;
ms = 80;

% Parâmetro AUC
u = 4;

analit_gammaBar_c = ones( length(analit_gammaBar) , max(N));
analit_gammaBar_c(:, 1) = analit_gammaBar; % variar só do primeiro canal...
simu_gammaBar_c = ones( length(simu_gammaBar) , max(N));
simu_gammaBar_c(:, 1) = simu_gammaBar; % variar só do primeiro canal...

% Erro de apontamento
% z = [0.6, 1.0, 1.5, 13];
% z = [0.6, 0.8; 1.0, 1.1; 3, 3.3; 14, 15];
% z = [0.6, 0.8; 1.0, 1.1; 1.5, 1.6; 3, 3.3; 14, 15];
% z = [0.6, 0.8; 1.0, 1.1; 1.5, 1.6; 14, 15];
z = [15];

for i=1:length(alpha)
    if ms <= (4/alpha(i))
        fprintf('ms > 4/alpha not met\n')
        error
    end
end

% Parâmetros da simulação
Nc = 1e3; % Número de pontos
rc = 1;   % 

% Perda de percurso
Hl = 1.00;

% Inicialização dos vetores -- Prealocation
% Pout = zeros(length(gammaBar));
% Ao = zeros(length(gammaBar));

colorz = 'brgmp';

sim_AUC = ones(1, length(simu_gammaBar));
sim_AUC = simulation_AUC(Nc, simu_gammaBar);

for w = 1:length(z)
    analit_params = [alpha, mu, ms, z(w);];
    
    [auc] = AUC_analit(N, analit_params, u, analit_gammaBar_c);
    % [gammaBar_dB, auc_asymp] = AUC_asymptotic(N, alpha, mu, ms, bounds, points, z(w, :), u);

    figure(1)
    % plot(gammaBar_dB, auc, colorz(w),...
    %      gammaBar_dB, auc_asymp,'k--',...
    %      'linewidth', 1.2)

    plot(analit_gammaBar_dB, auc, colorz(w),...
         simu_gammaBar_dB, sim_AUC, 'ok',...  
         'linewidth', 1.2)    
    
    hold on
end

axis([min(simu_gammaBar_dB) max(simu_gammaBar_dB) 0.5 1])
tam_fonte = 11;
legend();
% legend('FontSize', tam_fonte)
% legend("$z_1="+num2str(z(1,1))+", z_2="+num2str(z(1,2))+"$", '', "$z_1="+num2str(z(2,1))+", z_2="+num2str(z(2,2))+"$", '', "$z_1="+num2str(z(3,1))+", z_2="+num2str(z(3,2))+"$", '', "$z_1="+num2str(z(4,1))+", z_2="+num2str(z(4,2))+"$", 'Asymptotic', 'Location', 'southwest')
% legend("$z_1="+num2str(z(1,1))+", z_2="+num2str(z(1,2))+"$", '',"$z_1="+num2str(z(2,1))+", z_2="+num2str(z(2,2))+"$", '',"$z_1="+num2str(z(3,1))+", z_2="+num2str(z(3,2))+"$", '', "Non-pointing errors", 'Asymptotic', 'Location', 'southeast')
ax = gca;
ax.FontSize = 13;

% %title('Outage Probability')
% set(legend, 'Interpreter', 'latex')
% ylabel("Average AUC", 'FontSize', 14)
% % yticks([1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1e0])
% xlabel("SNR (dB)", 'FontSize', 14)
grid on

% %textbox com valores
% dim = [0.75 0.41 0.2 0.2];
% str = {"$\alpha_1 ="+num2str(alpha(1))+"$", "$\alpha_2 ="+num2str(alpha(2))+"$","$\mu ="+num2str(mu)+"$", "$m ="+num2str(ms)+"$", "$u ="+num2str(u)+"$"};
% annotation('textbox',dim,'interpreter','latex','String',str,'FitBoxToText','on', 'FontSize', tam_fonte);