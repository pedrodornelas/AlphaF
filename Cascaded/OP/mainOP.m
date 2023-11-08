
clear all
close all
clc

% Piece of code to calculate the Outage Probability
% Pedro Henrique Dornelas Almeida - 05/05/2023

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
simu_gammaBar_dB = linspace(L, U, 15); % SNR em dB
simu_gammaBar = 10.^(0.1*simu_gammaBar_dB); % SNR linear

% Limiares
% GthdB = linspace(-5,10,2);
% gth = 10.^(0.1*GthdB); % Limiar de SNR
gamma_th_dB = 5; %in dB
gamma_th = db2pow(gamma_th_dB);

% Parâmetros da Distribuição Alfa F - Cascaded
% N = [1, 2];    % nº estações relay
N = [1,2];    % nº estações relay
alfa = 2.5;
mu = [1.5, 1.7, 2];
ms = 2;

analit_gammaBar_c = 10.^(0.1 * ones( length(analit_gammaBar) , max(N)));
analit_gammaBar_c(:, 1) = analit_gammaBar; % variar só do primeiro canal...
simu_gammaBar_c = 10.^(0.1 * ones( length(simu_gammaBar) , max(N)));
simu_gammaBar_c(:, 1) = simu_gammaBar; % variar só do primeiro canal...

% Erro de apontamento
z = [0.6, 0.8;
     1  , 1.1;
     1.5, 1.6;
     8 , 9 ;];

% z = [0.6, 0.8;
%      1  , 1.1;
%      1.5, 1.6;
%      15 , 16 ;];
% z = [1, 1.1, 1.2];
% z = [15.0, 15.5, 16.0, 16.5];

if ms <= (4/alfa)
    fprintf('ms > 4/alpha not met')
    error
end

% Parâmetros da simulação
Nc = 1e6; % Número de pontos
rc = 1;   % hat r

% Perda de percurso
Hl = 1.00;

% Inicialização dos vetores -- Prealocation
Pout = zeros(1, length(simu_gammaBar));
% Ao = zeros(length(simu_gammaBar));

colorx = 'rb';
colorz = 'brgmp';
styleN = 'xo';

% for j = 1:length(N)
%     for k = 1:length(z)

tic;

figure(1)
h = [];
cont = 1;

for k = 1:size(z, 1)
    simulation_params = [alfa, mu(1), ms, z(k, 1), rc, Hl;
                         alfa, mu(2), ms, z(k, 2), rc, Hl;];

    analit_params = [alfa, mu(1), ms, z(k, 1);
                     alfa, mu(2), ms, z(k, 2);];

    gain_channels = individual_gain(max(N), simulation_params, Nc, simu_gammaBar_c);
    gain_total = cascaded_gain(gain_channels);
    
    % total_gain(:,:,2) = gain_channels(:,:,1) .* gain_channels(:,:,2) .* gain_channels(:,:,3)

    for j = 1:length(N)
        % [k j]
        % Gain = general_gain_cascaded(N(j), simulation_params, Nc, simu_gammaBar);
        Gain = gain_total(:,:,N(j)).';
        % size(Gain)
        for i = 1:length(simu_gammaBar)
            flagOP = sum(Gain(i,:).^2 <= gamma_th);
            Pout(i) = flagOP/Nc;
        end
        % Pout = ones(1, length(simu_gammaBar_dB));

        OP = OP_analit(N(j), analit_params, gamma_th, analit_gammaBar_c);
        OP_asy = OP_asymptotic(N(j), analit_params, gamma_th, analit_gammaBar_c);
        
        h(cont)   = semilogy(simu_gammaBar_dB, Pout, 'rx', 'linewidth', 1.2);hold on; % [colorz(j),'x']
        h(cont+1) = semilogy(analit_gammaBar_dB, OP, colorz(k), 'linewidth', 1.2);hold on;
        h(cont+2) = semilogy(analit_gammaBar_dB, OP_asy,'k--','linewidth', 1.2);hold on;

        cont = cont + 3;

        % semilogy(analit_gammaBar_dB, OP, [colorz(k)],...
        %          analit_gammaBar_dB, OP_asy,'k--',...
        %          'linewidth', 1.2)
        % hold on
    end
end

execution_time = toc;
disp(['Execution time: ' num2str(execution_time) 's']);


legend('FontSize', 11, 'Location', 'southwest')
% legend('' , "$z_1="+num2str(z(1,1))+", z_2="+num2str(z(1,2))+"$", '' , '' , '' , '' , '' , "$z_1="+num2str(z(2,1))+", z_2="+num2str(z(2,2))+"$", '' , '' , '' , '' , '' , "$z_1="+num2str(z(3,1))+", z_2="+num2str(z(3,2))+"$", '' , '' , '' , '' , '' ,  "Non-pointing errors",'Asymptotic', 'Simulated', 'Location', 'southwest')
legend([h(2), h(8), h(14), h(20), h(21), h(22)], { "$z_1="+num2str(z(1,1))+", z_2="+num2str(z(1,2))+"$" , "$z_1="+num2str(z(2,1))+", z_2="+num2str(z(2,2))+"$" , ...
                                                   "$z_1="+num2str(z(3,1))+", z_2="+num2str(z(3,2))+"$" , "Non-pointing errors", "Asymptotic", "Simulated" })

ax = gca;
ax.FontSize = 13;
tam_fonte=11;

axis([min(simu_gammaBar_dB) max(simu_gammaBar_dB) 1e-5 1])
set(legend, 'Interpreter', 'latex', 'LimitMaxLegendEntries', false)
ylabel("OP", 'FontSize', 14)
yticks([1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1e0])
xlabel("SNR", 'FontSize', 14)
grid on

%textbox com valores
dim = [0.15 0.45 0.2 0.2];
str = {"$\alpha ="+num2str(alfa)+"$", "$\mu_1 ="+num2str(mu(1))+"$", "$\mu_2 ="+num2str(mu(2))+"$" , "$m ="+num2str(ms)+"$", "$\gamma_{\rm th}="+num2str(gamma_th_dB)+"$ dB"};
annotation('textbox',dim,'interpreter','latex','String',str,'FitBoxToText','on', 'FontSize', tam_fonte);

%text arrow
% X = [0.74,0.7];
X = [0.62, 0.7];
% Y = [0.84,0.71];
Y = [0.19, 0.34];
str = {"$N=\{" + num2str(N(1)) + "," + num2str(N(2)) + "\}$"};
annotation('textarrow', X, Y, 'String', str, 'interpreter', 'latex', 'FontSize', 13);
