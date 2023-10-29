clear all
close all
clc

% Piece of code to calculate the Outage Probability
% Pedro Henrique Dornelas Almeida - 05/05/2023

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
N = [1, 2];    % nº estações relay
alfa = 2.5;
mu = [1.5, 1.7, 2];
ms = 2.0;

% Erro de apontamento
z = [1.0, 1.1];

if ms <= (4/alfa)
    fprintf('ms > 4/alpha not met')
    error
end

% Parâmetros da simulação
Nc = 1e4; % Número de pontos
rc = 1;   % hat r

% Perda de percurso
Hl = 1.00;

% Inicialização dos vetores -- Prealocation
Pout = zeros(1, length(simu_gammaBar));
% Ao = zeros(length(simu_gammaBar));

colorx = 'rb'
colorz = 'brgmp';

% for j = 1:length(N)
%     for k = 1:length(z)

simulation_params = [alfa, mu(1), ms, z(1), rc, Nc, Hl;
                     alfa, mu(2), ms, z(2), rc, Nc, Hl;];

analit_params = [alfa, mu(1), ms, z(1);
                 alfa, mu(2), ms, z(2);];

for j = 1:length(N)
    % params = [alfa, mu(1), ms, z(1), rc, Nc, Hl];
    Gain = general_gain_cascaded(N(j), simulation_params, simu_gammaBar);
    % line = size(Gain,1);
    % col = size(Gain,2);
    for i = 1:length(simu_gammaBar)
        flagOP = sum(Gain(i,:).^2 <= gamma_th);
        Pout(i) = flagOP/Nc;
    end

    OP = OP_analit(N(j), analit_params, gamma_th, analit_gammaBar);

    figure(1)
    semilogy(simu_gammaBar_dB, Pout, [colorz(j),'x'],...
             analit_gammaBar_dB, OP, colorz(j),...
             ...%  gammaBar_dB, P,'k--',...
             'linewidth',1.2)
    hold on
end

        % Ao = sqrt(simu_gammaBar*(2+z(k)^2))/(rc*z(k)*Hl);
        % for i = 1:length(simu_gammaBar)
        %     [j k i]
        %     % Ganhos aleatórios
        %     Hf = gainAF(alfa(j),mu,ms,rc,Nc,-1e-3); % Gain Alpha F
        %     Hp = PointError(z(k),Ao(i),Nc);         % Gain Pointing error
        %     % Ganho total
        %     Gain = (Hl(:).*Hf(:).*Hp(:)).';
            
        %     [flagOP,~] = find(Gain.^2 <= gamma_th);
        %     Pout(i) = sum(flagOP)/Nc;
        % end

        %%
        
        % nj = N(j);
        % [gammaBar_dB, Pb] = OP_analit(N(j), alfa, mu(1:nj), ms, bounds, points, z(k, 1:nj), gamma_th);
        % [gammaBar_dB, P] = OP_asymptotic(N(j), alfa, mu(1:nj), ms, bounds, points, z(k, 1:nj), gamma_th);

        % figure(1)
        % semilogy(simu_gammaBar_dB, Pout(:,1),'rx',...
        %         gammaBar_dB, Pb, colorz(j),...
        %         gammaBar_dB, P,'k--',...
        %         'linewidth',1.2)
        % str = colorz(j)+'o'
        
        % if k == 4
        %     color = 'm--'
        % else
        %     color = colorz(k)
        % end
        % semilogy(gammaBar_dB, Pb, colorz(k),...
        %          gammaBar_dB, P,'k--',...
        %          'linewidth',1.2)
        
        % hold on
%     end
% end

axis([min(simu_gammaBar_dB) max(simu_gammaBar_dB) 1e-5 1])
set(legend, 'Interpreter', 'latex')
ylabel("OP", 'FontSize', 14)
yticks([1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1e0])
xlabel("SNR", 'FontSize', 14)
grid on
