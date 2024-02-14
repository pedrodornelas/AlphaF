clc
clear all;
% close all;

% A piece of code to calculate the Bit Error Probability (BEP)
% Pedro Henrique Dornelas Almeida - 14/02/2023

functions_path = "functions";
addpath(functions_path);

set(0,'defaulttextinterpreter','latex');
set(groot,'defaultAxesTickLabelInterpreter','latex');

load('BEP.mat');

% distribution params
L = [1,2,3,4];
N = 2;
alpha = [1.5, 2.3];

ms = [4, 5];
mu = [2.5, 3.5];
z = [0.7, 0.8;
       7, 7.2;];

p = 1;  % modulation param (p=1: BPSK)
M = 2;  % order constelation MQAM

% simulation params
rc = 1;
Hl = 1;
Nc = 5e6;

% SNR
L_bound_dB = 0;
U_bound_dB = 30;
points = 15;
simu_gammaBar_dB = linspace(L_bound_dB, U_bound_dB, points);
simu_gammaBar = 10 .^ (0.1 * simu_gammaBar_dB);
% individual channels gammaBar
simu_gammaBar_c = db2pow(1) * ones(length(simu_gammaBar), max(N));
simu_gammaBar_c(:, 1) = simu_gammaBar; % variar somente primeiro canal

tic;
figure()
colorz = 'brgmp';
cont = 1;
h = [];

for j = 1:length(z)
    simulation_params = [alpha(1), mu(1), ms(1), z(j, 1), rc, Hl;
                         alpha(2), mu(2), ms(2), z(j, 2), rc, Hl;];

    BEP_simulated = channel(M, N, L, simulation_params, Nc, simu_gammaBar_c);

    for i = 1:length(L)
        h(cont)   = semilogy(simu_gammaBar_dB, BEP_simulated(:, i), 'rx', 'linewidth', 1.2); hold on;
        h(cont+1) = semilogy(gamma_bar_dB, BEP(:, i, j), colorz(i), 'linewidth', 1.2); hold on;
        h(cont+2) = semilogy(gamma_bar_dB, BEP_asy(:, i, j), 'k--', 'linewidth', 1.2); hold on;
        cont = cont+3;
    end
end

execution_time = toc;
format_time = executionTime(execution_time);
disp(['Execution time: ' format_time]);

axis([min(gamma_bar_dB) max(gamma_bar_dB) 1e-5 1])
tam_fonte = 10;
legend('FontSize', tam_fonte, 'Location', 'southeast')
legend([h(2), h(5), h(8), h(11), h(12), h(13)], {"$L="+num2str(L(1))+"$", "$L="+num2str(L(2))+"$", "$L="+num2str(L(3))+"$", "$L="+num2str(L(4))+"$", "Asymptotic", "Simulated"});
% legend([h(1), h(5), h(9), h(10)], {"$L="+num2str(L(1))+"$", "$L="+num2str(L(2))+"$", "$L="+num2str(L(3))+"$", "Asymptotic"});
ax = gca;
ax.FontSize = 13;

%title('Outage Probability')
set(legend, 'Interpreter', 'latex')
ylabel("BEP", 'FontSize', 14)
yticks([1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1e0])
xlabel("SNR (dB)", 'FontSize', 14)
grid on