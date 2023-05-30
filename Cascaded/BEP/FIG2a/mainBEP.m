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

GBdB = linspace(L, U, 15); % SNR em dB
gammaBar = 10.^(0.1*GBdB); % SNR linear

% Parâmetro da modulação
rho = 1;

% Parâmetros da Distribuição Alfa F - Cascaded
N = 2;    % nº estações relay
alfa = [2.2, 2.5];
mu = 1;
ms = [2];

% Erro de apontamento
z = [0.6, 0.8; 1.0, 1.1; 1.5, 1.6; 14, 15];
% z = [0.6, 0.8; 1.0, 1.1; 3, 3.3; 14, 15];

for i=1:length(alfa)
    if ms <= (4/alfa(i))
        fprintf('ms > 4/alpha not met\n')
        error
    end
end

% Parâmetros da simulação
Nc = 1e4; % Número de pontos
rc = 1;   % 

% Perda de percurso
Hl = 1.00;

% Inicialização dos vetores -- Prealocation
Pout = zeros(length(gammaBar));
Ao = zeros(length(gammaBar));

colorz = 'brgmp';

for j = 1:length(alfa)
    for k = 1:length(z)
        for w = 1:length(ms)
            % Ao = sqrt(gammaBar*(2+z(k)^2))/(rc*z(k)*Hl);
            % for i = 1:length(gammaBar)
            %     [j k i]
            %     % Ganhos aleatórios
            %     Hf = gainAF(alfa(j),mu,ms,rc,Nc,-1e-3); % Alpha F
            %     Hp = PointError(z(k),Ao(i),Nc); % Pointing error
            %     % Ganho total
            %     Gain = (Hl(:).*Hf(:).*Hp(:)).';
                
            %     [flagBEP,~] = find(Gain.^2 <= rho);
            %     Pout(i) = sum(flagBEP)/Nc;
            % end

            %%
            
            [gammaBar_dB, Pb] = BEP_analit(N, alfa(1:N), mu, ms(w), bounds, points, z(k, 1:N), rho);
            [gammaBar_dB, P] = BEP_asymptotic(N, alfa(1:N), mu, ms(w), bounds, points, z(k, 1:N), rho);

            figure(1)
            % semilogy(GBdB, Pout(:,1),'rx',...
            %         gammaBar_dB, Pb, colorz(j),...
            %         gammaBar_dB, P,'k--',...
            %         'linewidth',1.2)
            % str = colorz(j)+'o'
            
            % if k == 4
            %     color = 'm--'
            % else
            %     color = colorz(k)
            % end
            semilogy(gammaBar_dB, Pb, colorz(k),...
                    gammaBar_dB, P,'k--',...
                    'linewidth',1.2)
            
            hold on
        end
    end
end

axis([min(GBdB) max(GBdB) 1e-5 1])
tam_fonte = 11;
legend('FontSize', tam_fonte)
% legend("$z_1="+num2str(z(1,1))+", z_2="+num2str(z(1,2))+"$", '', "$z_1="+num2str(z(2,1))+", z_2="+num2str(z(2,2))+"$", '', "$z_1="+num2str(z(3,1))+", z_2="+num2str(z(3,2))+"$", '', "$z_1="+num2str(z(4,1))+", z_2="+num2str(z(4,2))+"$", 'Asymptotic', 'Location', 'southwest')
legend("$z_1="+num2str(z(1,1))+", z_2="+num2str(z(1,2))+"$", '', "$z_1="+num2str(z(2,1))+", z_2="+num2str(z(2,2))+"$", '', "$z_1="+num2str(z(3,1))+", z_2="+num2str(z(3,2))+"$", '', "Non-pointing errors", 'Asymptotic', 'Location', 'southwest')
ax = gca;
ax.FontSize = 13;

%title('Outage Probability')
set(legend, 'Interpreter', 'latex')
ylabel("BEP", 'FontSize', 14)
yticks([1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1e0])
xlabel("SNR (dB)", 'FontSize', 14)
grid on

if rho == 1
    modulation = 'BPSK';
elseif rho == 0.5
    modulation = 'QPSK';
end

%textbox com valores
dim = [0.15 0.32 0.2 0.2];
str = {"$\alpha_1 ="+num2str(alfa(1))+"$", "$\alpha_2 ="+num2str(alfa(2))+"$", "$\mu ="+num2str(mu)+"$", "$m ="+num2str(ms)+"$","$\rho="+num2str(rho)+"$ (BSPK)"};
annotation('textbox',dim,'interpreter','latex','String',str,'FitBoxToText','on', 'FontSize', tam_fonte);

%GBdB,cdfG(gth(1),gammaBar),'-',