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

% Parâmetros da Distribuição Alpha F - Cascaded
N = 2;              % nº estações relay
alpha = [2.0, 3.7];
mu = 1;
ms = 2.3;

% Erro de apontamento
z = [0.6, 1.0, 1.5, 13];
% z = [0.6, 0.8; 1.0, 1.1; 3, 3.3; 14, 15];

for i=1:length(alpha)
    if ms <= (4/alpha(i))
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

for j = 1:length(alpha)
    for k = 1:length(z)
        % Ao = sqrt(gammaBar*(2+z(k)^2))/(rc*z(k)*Hl);
        % for i = 1:length(gammaBar)
        %     [j k i]
        %     % Ganhos aleatórios
        %     Hf = gainAF(alpha(j),mu,ms,rc,Nc,-1e-3); % Alpha F
        %     Hp = PointError(z(k),Ao(i),Nc); % Pointing error
        %     % Ganho total
        %     Gain = (Hl(:).*Hf(:).*Hp(:)).';
                
        %     [flagCAP,~] = find(Gain.^2 <= rho);
        %     Pout(i) = sum(flagCAP)/Nc;
        % end

        %%
            
        [gammaBar_dB, Pb] = CAP_analit(N, alpha(j), mu, ms, bounds, points, z(k));
        [gammaBar_dB, P] = CAP_asymptotic(N, alpha(j), mu, ms, bounds, points, z(k));

        figure(1)
        % plot(GBdB, Pout(:,1),'rx',...
        %         gammaBar_dB, Pb, colorz(j),...
        %         gammaBar_dB, P,'k--',...
        %         'linewidth',1.2)
        % str = colorz(j)+'o'
            
        % if k == 4
        %     color = 'm--'
        % else
        %     color = colorz(k)
        % end
        plot(gammaBar_dB, Pb, colorz(k),...
             gammaBar_dB, P,'k--',...
             'linewidth', 1.2)
        
        hold on
    end
end

axis([min(GBdB) max(GBdB) 0 20])
tam_fonte = 11;
legend('FontSize', tam_fonte)
% legend("$z_1="+num2str(z(1,1))+", z_2="+num2str(z(1,2))+"$", '', "$z_1="+num2str(z(2,1))+", z_2="+num2str(z(2,2))+"$", '', "$z_1="+num2str(z(3,1))+", z_2="+num2str(z(3,2))+"$", '', "$z_1="+num2str(z(4,1))+", z_2="+num2str(z(4,2))+"$", 'Asymptotic', 'Location', 'southwest')
legend("$z="+num2str(z(1))+"$", '', "$z="+num2str(z(2))+"$", '', "$z="+num2str(z(3))+"$", '', "Non-pointing errors", 'Asymptotic', 'Location', 'northwest')
ax = gca;
ax.FontSize = 13;

%title('Outage Probability')
set(legend, 'Interpreter', 'latex')
ylabel("Capacity", 'FontSize', 14)
% yticks([1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1e0])
xlabel("SNR (dB)", 'FontSize', 14)
grid on

%textbox com valores
dim = [0.15 0.45 0.2 0.2];
str = {"$\mu ="+num2str(mu)+"$", "$m ="+num2str(ms)+"$"};
annotation('textbox',dim,'interpreter','latex','String',str,'FitBoxToText','on', 'FontSize', tam_fonte);

%arrow
X = [0.69,0.76];
Y = [0.74,0.56];
annotation('arrow', X, Y);

% text
X = 37.47;
Y = 10.21;
str = {"$\alpha=\{"+num2str(alpha(2))+","+num2str(alpha(1))+"\}$"};
text(X, Y, str, 'Interpreter', 'latex','FontSize', tam_fonte);

% 35.54, 4.20, 1.4e
% 35.5466,3.8197,7.9179,0.7801
%35.63868103130755,4.265417730496456,1.4e
