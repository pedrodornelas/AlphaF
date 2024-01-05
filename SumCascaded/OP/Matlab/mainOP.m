clear all
close all

% A piece of code to calculate the Outage Probability
% Pedro Henrique Dornelas Almeida - 03/11/2023

functions_path = "functions";
addpath(functions_path);

verify_python();
% Se o caminho do Python foi encontrado, continue com o código
disp('Continuing with code...');


set(0,'defaulttextinterpreter','latex');
set(groot,'defaultAxesTickLabelInterpreter','latex');


% distribution parameters
% L = 1; % sum of cascaded channels
% N = [1,2]; % cascaded channels
% ms = [2];
% alpha = [2.5];
% mu = [1.5, 1.7];
% z = [0.7, 0.8;
%      1  , 1.1;
%      1.5, 1.6;
%      8 , 9 ;];

% RIS a-F with pointing errors
% L = [1,2,3,4]; % sum of cascaded channels
L = [1,2]
N = [2]; % cascaded channels
ms = [3,4];
alpha = [1.5, 2.3];
mu = [1.5, 1.7];
z = [0.7, 0.8;
     7  , 8;];

% z = [7, 8;];

% simulation parameters
rc = 1;
Hl = 1;
Nc = 1e5;


% threshold OP
gamma_th_dB = 5; %in dB
gamma_th = db2pow(gamma_th_dB);

for j = 1:length(alpha)
    % j
    if ms <= (4/alpha(j))
        error("ms > 4/alpha not met. Exiting...")
    end
end

% SNRs -- Amostragem dos valores observáveis
L_bound = 0;    %db
U_bound = 50;   %db
points = 1e2;
bounds = [L_bound U_bound]; %dB gammaBar limits

analit_gammaBar_dB = linspace(L_bound, U_bound, points);
analit_gammaBar = 10.^(0.1*analit_gammaBar_dB);
simu_gammaBar_dB = linspace(L_bound, U_bound, 15); % SNR em dB
simu_gammaBar = 10.^(0.1*simu_gammaBar_dB); % SNR linear

% individual channels gammaBar
analit_gammaBar_c = db2pow(1) * ones( length(analit_gammaBar) , max(N));
analit_gammaBar_c(:, 1) = analit_gammaBar; % variar só do primeiro canal...
simu_gammaBar_c = db2pow(1) * ones( length(simu_gammaBar) , max(N));
simu_gammaBar_c(:, 1) = simu_gammaBar; % variar só do primeiro canal...

% pre-allocatin vector OP
Pout = zeros(1, length(simu_gammaBar));

tic;
colorz = 'brgmp';
stylez = 'xo';
h = [];
cont = 1;
figure()

for i = 1:length(L)
    % L(i)
    for k = 1:size(z, 1)
        % k
        analit_params = [alpha(1), mu(1), ms(1), z(k, 1);
                         alpha(2), mu(2), ms(2), z(k, 2);];

        simulation_params = [alpha(1), mu(1), ms(1), z(k, 1), rc, Hl;
                             alpha(2), mu(2), ms(2), z(k, 2), rc, Hl;];
        
        % gain_channels = individual_gain(max(N), simulation_params, Nc, simu_gammaBar_c);
        % gain_cascaded = cascaded_gain(gain_channels);
        % gain_sum_cascaded = sum_cascaded_gain(L(i), N, gain_cascaded);

        for j = 1:length(N)
            [k, L(i)]

            % Gain = gain_sum_cascaded(:, :, L(i)).';

            % for k = 1:length(simu_gammaBar)
            %     flagOP = sum(Gain(k, :).^2 <= gamma_th);
            %     Pout(k) = flagOP/Nc;
            % end

            OP = OP_analit(L(i), N(j), analit_params, gamma_th, analit_gammaBar_c);
            OP_asy = OP_asymptotic(L(i), N(j), analit_params, gamma_th, analit_gammaBar_c);

            % h(cont) = semilogy(simu_gammaBar_dB, Pout, 'rx', 'linewidth', 1.2);hold on;
            % cont = cont + 1;
            h(cont) = semilogy(analit_gammaBar_dB, OP, colorz(i), 'linewidth', 1.2);hold on;
            cont = cont + 1;
            h(cont) = semilogy(analit_gammaBar_dB, OP_asy, ['k--'], 'linewidth', 1.2);hold on;
            cont = cont + 1;
        end
    end
end

execution_time = toc;
disp(['Execution time: ' num2str(execution_time) 's']);

ax = gca;
ax.FontSize = 13;
tam_fonte=11;

legend('FontSize', 11, 'Location', 'southwest')
% sum cascaded legend
% legend([h(1), h(3), h(5), h(7)], {"$z_1="+num2str(z(1,1))+", z_2="+num2str(z(1,2))+"$" , "$z_1="+num2str(z(2,1))+", z_2="+num2str(z(2,2))+"$" , ...
%                                   "$z_1="+num2str(z(3,1))+", z_2="+num2str(z(3,2))+"$" , "Non-pointing errors"})
% RIS legend
% legend([h(1), h(5), h(9), h(13), h(14)], {"$L="+num2str(L(1))+"$", "$L="+num2str(L(2))+"$", "$L="+num2str(L(3))+"$", "$L="+num2str(L(4))+"$", "Asymptotic"})
% legend([h(2), h(8), h(14), h(20), h(21), h(22)], {"$L="+num2str(L(1))+"$", "$L="+num2str(L(2))+"$", "$L="+num2str(L(3))+"$", "$L="+num2str(L(4))+"$", "Asymptotic", "Simulated"});
legend([h(1), h(5), h(6)], {"$L="+num2str(L(1))+"$", "$L="+num2str(L(2))+"$", "Asymptotic"})
% legend([h(1), h(3), h(5), h(7), h(8)], {"$L="+num2str(L(1))+"$", "$L="+num2str(L(2))+"$", "$L="+num2str(L(3))+"$", "$L="+num2str(L(4))+"$", "Asymptotic"})



axis([min(analit_gammaBar_dB) max(analit_gammaBar_dB) 1e-5 1])
set(legend, 'Interpreter', 'latex')
ylabel("OP", 'FontSize', 14)
yticks([1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1e0])
xlabel("SNR", 'FontSize', 14)
grid on

% textbox com valores
dim = [0.60,0.14,0.23,0.23];
% dim = [0.15 0.45 0.2 0.2];
str = {"$N="+num2str(N)+"$" , "$\alpha_1 ="+num2str(alpha(1))+", \alpha_2 ="+num2str(alpha(2))+"$", "$\mu_1 ="+num2str(mu(1))+", \mu_2 ="+num2str(mu(2))+"$" , "$m_1 ="+num2str(ms(1))+", m_2 ="+num2str(ms(2))+"$", "$\gamma_{\rm th}="+num2str(gamma_th_dB)+"$ dB"};
annotation('textbox',dim,'interpreter','latex','String',str,'FitBoxToText','on', 'FontSize', tam_fonte);
