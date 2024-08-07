clear all
close all
clc

% Piece of code to calculate the Capacity
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
simu_gammaBar_dB = linspace(L, U, 15); % SNR em dB
simu_gammaBar = 10.^(0.1*simu_gammaBar_dB); % SNR linear

% Parâmetros da Distribuição Alpha F with pointing error- Cascaded
N = 2;               % nº estações relay
alpha = [2.0, 3.7];  % non-linearity
mu = 1.8;            % multipath cluster
ms = 2.5;            % shadowing
% Erro de apontamento
z = [0.7, 1.0, 1.3, 8];
% z = [0.6, 0.8; 1.0, 1.1; 3, 3.3; 14, 15];

for i=1:length(alpha)
    if ms <= (4/alpha(i))
        fprintf('ms > 4/alpha not met\n')
        error
    end
end

% Parâmetros da simulação
Nc = 1e5;  % points number
rc = 1;    % hat r
Hl = 1.00; % Perda de percurso

% Inicialização dos vetores -- Prealocation
CAP_sim = zeros(1, length(simu_gammaBar));

analit_gammaBar_c = db2pow(1) * (ones( length(analit_gammaBar) , max(N)));
analit_gammaBar_c(:, 1) = analit_gammaBar; % variar só do primeiro canal...
simu_gammaBar_c = db2pow(1) * (ones( length(simu_gammaBar) , max(N)));
simu_gammaBar_c(:, 1) = simu_gammaBar; % variar só do primeiro canal...

colorz = 'brgmp';
styleA = 'xo';

tic;

% close(figure(1));
figure(1)
h = [];
cont = 1;


for j = 1:length(alpha)
    for k = 1:length(z)

        simulation_params = [alpha(j), mu, ms, z(k), rc, Hl;
                             alpha(j), mu, ms, z(k), rc, Hl;];
        
        analit_params = [alpha(j), mu, ms, z(k);
                         alpha(j), mu, ms, z(k);];

        gain_channels = individual_gain(max(N), simulation_params, Nc, simu_gammaBar_c);
        gain_total = cascaded_gain(gain_channels);

        Gain = gain_total(:,:, N);
        for i = 1:length(simu_gammaBar)
            CAP_sim(i) = mean(log2(1+Gain(:,i).^2)).';
        end
            
        CAP = CAP_analit(N, analit_params, analit_gammaBar_c);
        CAP_asy = CAP_asymptotic(N, analit_params, analit_gammaBar_c);

        h(cont)   = plot(analit_gammaBar_dB, CAP, [colorz(k)], 'linewidth', 1.2);hold on;
        h(cont+1) = plot(analit_gammaBar_dB, CAP_asy, 'k--', 'linewidth', 1.2);hold on;
        h(cont+2) = plot(simu_gammaBar_dB, CAP_sim, 'rx', 'linewidth', 1.2);hold on;

        cont = cont + 3;
    end
end

execution_time = toc;
format_time = executionTime(execution_time);
disp(['Execution time: ' format_time]);

axis([min(simu_gammaBar_dB) max(simu_gammaBar_dB) 0 10])
tam_fonte = 11;
legend('FontSize', tam_fonte, 'Location', 'northwest')
% legend("$z_1="+num2str(z(1,1))+", z_2="+num2str(z(1,2))+"$", '', "$z_1="+num2str(z(2,1))+", z_2="+num2str(z(2,2))+"$", '', "$z_1="+num2str(z(3,1))+", z_2="+num2str(z(3,2))+"$", '', "$z_1="+num2str(z(4,1))+", z_2="+num2str(z(4,2))+"$", 'Asymptotic', 'Location', 'southwest')
% legend("$z="+num2str(z(1))+"$", '', "$z="+num2str(z(2))+"$", '', "$z="+num2str(z(3))+"$", '', "Non-pointing errors", 'Asymptotic', 'Location', 'northwest')
legend([h(1), h(4), h(7), h(10), h(20), h(21)], {"$z="+num2str(z(1))+"$", "$z="+num2str(z(2))+"$", "$z="+num2str(z(3))+"$", ...
                                                  "Non-pointing errors", "Asymptotic", "Simulated"})
ax = gca;
ax.FontSize = 13;

%title('Outage Probability')
set(legend, 'Interpreter', 'latex')
ylabel("Capacity", 'FontSize', 14)
% yticks([1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1e0])
xlabel("SNR (dB)", 'FontSize', 14)
grid on

%textbox com valores
dim = [0.15 0.40 0.2 0.2];
str = {"$\mu ="+num2str(mu)+"$", "$m ="+num2str(ms)+"$"};
annotation('textbox',dim,'interpreter','latex','String',str,'FitBoxToText','on', 'FontSize', tam_fonte);

%arrow
X = [0.541071428571429,0.564285714285712];
Y = [0.452380952380953,0.345238095238099];
annotation('arrow', X, Y);

% text
X = 23.77198156682029;
Y = 2.415130125150542;
str = {"$\alpha=\{"+num2str(alpha(2))+","+num2str(alpha(1))+"\}$"};
text(X, Y, str, 'Interpreter', 'latex','FontSize', tam_fonte);

% 35.54, 4.20, 1.4e
% 35.5466,3.8197,7.9179,0.7801
%35.63868103130755,4.265417730496456,1.4e
