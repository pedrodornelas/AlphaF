clear all
close all
clc

% Piece of code to calculate the BEP - Bit Error Probability
% Pedro Henrique Dornelas Almeida - 08/11/2023

functions_path = "functions";
addpath(functions_path);

set(0,'defaulttextinterpreter','latex');
set(groot,'defaultAxesTickLabelInterpreter','latex');

% SNRs -- Amostragem dos valores observáveis
L = 0;    %db
U = 50;   %db
points = 1e2;
bounds = [L U]; %dB gammaBar limits

analit_gammaBar_dB = linspace(L, U, points);
analit_gammaBar = 10.^(0.1*analit_gammaBar_dB);
simu_gammaBar_dB = linspace(L, U, 15);      % SNR em dB
simu_gammaBar = 10.^(0.1*simu_gammaBar_dB); % SNR linear

% Parâmetro da modulação (teóricas)
rho = 1;
% Parâmetro da modulação (simulação)
% Ordem da constelação MQAM
M = 2;

% Parâmetros da Distribuição Alfa F - Cascaded
N = 2;    % nº estações relay
alfa = [2.2, 2.5];
mu = 1;
ms = [2];

analit_gammaBar_c = db2pow(1) * (ones( length(analit_gammaBar) , max(N)));
% analit_gammaBar_c = (ones( length(analit_gammaBar) , max(N)));
analit_gammaBar_c(:, 1) = analit_gammaBar; % variar só do primeiro canal...
simu_gammaBar_c = db2pow(1) * (ones( length(simu_gammaBar) , max(N)));
% simu_gammaBar_c = (ones( length(simu_gammaBar) , max(N)));
simu_gammaBar_c(:, 1) = simu_gammaBar; % variar só do primeiro canal...

% Erro de apontamento
z = [0.6, 0.8; 
     1.0, 1.1; 
     1.5, 1.6; 
     8 ,  9;];
% z = [ 8,  9;
%      14, 15;];
% z = [0.6, 0.8; 1.0, 1.1; 3, 3.3; 14, 15];

for i=1:length(alfa)
    if ms <= (4/alfa(i))
        fprintf('ms > 4/alpha not met\n')
        error
    end
end

% Parâmetros da simulação
Nc = 1e6; % Número de pontos
Ao = 0.8;
% Ao = 1;
rc = 1;

% Perda de percurso
Hl = 1.00;

% início da simulação e geração de gráficos
colorz = 'brgmp';

tic;


figure(1)
cont = 1;


for k = 1:length(z)
    k
    for w = 1:length(ms)
        % simulation_params = [alfa(1), mu, ms(w), z(k, 1), rc, Hl;
        %                      alfa(2), mu, ms(w), z(k, 2), rc, Hl;];

        simulation_params = [alfa(1), mu, ms(w), z(k, 1), Ao, Hl;
                             alfa(2), mu, ms(w), z(k, 2), Ao, Hl;];

        analit_params = [alfa(1), mu, ms(w), z(k, 1);
                         alfa(2), mu, ms(w), z(k, 2);];

        BEP_sim = channel(M, N, simulation_params, Nc, simu_gammaBar_c);        
        [BEP] = BEP_analit(N, analit_params, rho, analit_gammaBar_c);
        [BEP_asy] = BEP_asymptotic(N, analit_params, rho, analit_gammaBar_c);

        h(cont)   = semilogy(analit_gammaBar_dB, BEP, colorz(k),'linewidth',1.2); hold on;
        h(cont+1) = semilogy(analit_gammaBar_dB, BEP_asy,'k--', 'linewidth',1.2);hold on;
        h(cont+2) = semilogy(simu_gammaBar_dB, BEP_sim, [colorz(k),'x'], 'linewidth', 1.2);hold on;
        cont = cont+3;
        
    end
end

execution_time = toc;
disp(['Execution time: ' num2str(execution_time) 's']);

axis([min(simu_gammaBar_dB) max(simu_gammaBar_dB) 1e-5 1])
tam_fonte = 11;
legend('FontSize', tam_fonte, 'Location', 'southwest')
% legend("$z_1="+num2str(z(1,1))+", z_2="+num2str(z(1,2))+"$", '', "$z_1="+num2str(z(2,1))+", z_2="+num2str(z(2,2))+"$", '', "$z_1="+num2str(z(3,1))+", z_2="+num2str(z(3,2))+"$", '', "$z_1="+num2str(z(4,1))+", z_2="+num2str(z(4,2))+"$", 'Asymptotic', 'Location', 'southwest')
% legend("$z_1="+num2str(z(1,1))+", z_2="+num2str(z(1,2))+"$", '', "$z_1="+num2str(z(2,1))+", z_2="+num2str(z(2,2))+"$", '', "$z_1="+num2str(z(3,1))+", z_2="+num2str(z(3,2))+"$", '', "Non-pointing errors", 'Asymptotic', 'Location', 'southwest')
legend([h(1), h(4), h(7), h(10), h(11), h(12)], {"$z_1="+num2str(z(1,1))+", z_2="+num2str(z(1,2))+"$" , "$z_1="+num2str(z(2,1))+", z_2="+num2str(z(2,2))+"$",...
                                                 "$z_1="+num2str(z(3,1))+", z_2="+num2str(z(3,2))+"$" , "Non-pointing errors", "Asymptotic", "Simulated"});
ax = gca;
ax.FontSize = 13;

%title('Outage Probability')
set(legend, 'Interpreter', 'latex')
ylabel("BEP", 'FontSize', 14)
yticks([1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1e0])
xlabel("SNR (dB)", 'FontSize', 14)
grid on

if rho == 1
    modulation = '(BPSK)';
elseif rho == 0.5
    modulation = '(QPSK)';
end

%textbox com valores
dim = [0.15 0.45 0.2 0.2];
str = {"$\alpha_1 ="+num2str(alfa(1))+"$", "$\alpha_2 ="+num2str(alfa(2))+"$", "$\mu ="+num2str(mu)+"$", "$m ="+num2str(ms)+"$","$\rho="+num2str(rho)+"$ "+modulation};
annotation('textbox',dim,'interpreter','latex','String',str,'FitBoxToText','on', 'FontSize', tam_fonte);

%simu_gammaBar_dB,cdfG(gth(1),gammaBar),'-',