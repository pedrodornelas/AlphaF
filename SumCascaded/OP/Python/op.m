clc
clear all;
% close all;

% A piece of code to calculate the Outage Probability
% Pedro Henrique Dornelas Almeida - 06/02/2023

functions_path = "functions";
addpath(functions_path);

set(0,'defaulttextinterpreter','latex');
set(groot,'defaultAxesTickLabelInterpreter','latex');

load('OP.mat');

L = [1,2,3,4];
N = 2;
alpha = [1.5, 2.3];
ms = [3, 4];
mu = [1.5, 1.7];
z = [0.7, 0.8;
       7,   8;];


% simulation params
rc = 1;
Hl = 1;
Nc = 5e5;
% Nc = 1e3;

% threshold OP
gamma_th_dB = 5;
gamma_th = db2pow(gamma_th_dB);

% SNR
L_bound_dB = 0;
U_bound_dB = 50;
points = 15;
simu_gammaBar_dB = linspace(L_bound_dB, U_bound_dB, points);
simu_gammaBar = 10 .^ (0.1 * simu_gammaBar_dB);
% individual channels gammaBar
simu_gammaBar_c = db2pow(1) * ones(length(simu_gammaBar), max(N));
simu_gammaBar_c(:, 1) = simu_gammaBar; % variar somente primeiro canal

Pout = zeros(1, length(simu_gammaBar));

tic;
figure()
colorz = 'brgmp';
cont = 1;
h = [];

for j = 1:length(z)
    simulation_params = [alpha(1), mu(1), ms(1), z(j, 1), rc, Hl;
                         alpha(2), mu(2), ms(2), z(j, 2), rc, Hl;];
        
    gain_sum_cascaded = sum_cascaded_gain(L, N, Nc, points, simu_gammaBar_c, simulation_params);
    for i = 1:(length(L))
        
        gain_sum = gain_sum_cascaded(:,:,i);
        for k = 1:length(simu_gammaBar)
            flagOP = sum(gain_sum(:, k).^2 <= gamma_th);
            Pout(k) = flagOP/Nc;
        end

        h(cont)   = semilogy(simu_gammaBar_dB, Pout, 'rx', 'linewidth', 1.2);hold on;
        h(cont+1) = semilogy(gamma_bar_dB, OP(:, i, j), colorz(i), 'linewidth', 1.2); hold on;
        h(cont+2) = semilogy(gamma_bar_dB, OP_asy(:, i, j), 'k--', 'linewidth', 1.2); hold on;
        cont = cont+3;
    end
end

execution_time = toc;
disp(['Execution time: ' num2str(execution_time) 's']);


axis([min(gamma_bar_dB) max(gamma_bar_dB) 1e-5 1])
tam_fonte = 11;
legend('FontSize', tam_fonte, 'Location', 'southwest')
legend([h(2), h(5), h(8), h(11), h(12), h(13)], {"$L="+num2str(L(1))+"$", "$L="+num2str(L(2))+"$", "$L="+num2str(L(3))+"$", "$L="+num2str(L(4))+"$", "Asymptotic", "Simulated"});
% legend([h(1), h(5), h(9), h(13), h(14)], {"$L="+num2str(L(1))+"$", "$L="+num2str(L(2))+"$", "$L="+num2str(L(3))+"$", "$L="+num2str(L(4))+"$", "Asymptotic"});
% legend([h(1), h(5), h(9), h(10)], {"$L="+num2str(L(1))+"$", "$L="+num2str(L(2))+"$", "$L="+num2str(L(3))+"$", "Asymptotic"});
ax = gca;
ax.FontSize = 13;

%title('Outage Probability')
set(legend, 'Interpreter', 'latex')
ylabel("OP", 'FontSize', 14)
yticks([1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1e0])
xlabel("SNR (dB)", 'FontSize', 14)
grid on

% textbox com valores
dim = [0.60,0.14,0.23,0.23];
% dim = [0.15 0.45 0.2 0.2];
str = {"$N="+num2str(N)+"$" , "$\alpha_1 ="+num2str(alpha(1))+", \alpha_2 ="+num2str(alpha(2))+"$", "$\mu_1 ="+num2str(mu(1))+", \mu_2 ="+num2str(mu(2))+"$" , "$m_1 ="+num2str(ms(1))+", m_2 ="+num2str(ms(2))+"$", "$\gamma_{\rm th}="+num2str(gamma_th_dB)+"$ dB"};
annotation('textbox',dim,'interpreter','latex','String',str,'FitBoxToText','on', 'FontSize', tam_fonte);