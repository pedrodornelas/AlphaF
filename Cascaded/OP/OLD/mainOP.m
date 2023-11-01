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

% Limiares
% GthdB = linspace(-5,10,2);
% gth = 10.^(0.1*GthdB); % Limiar de SNR
gamma_th_dB = 5; %in dB
gamma_th = db2pow(gamma_th_dB);

% Parâmetros da Distribuição Alfa F - Cascaded
N = [1, 2];    % nº estações relay
alfa = 2.5;
mu = [1.5, 1.7, 2];
ms = 2.0;

% Erro de apontamento
z = [0.6, 0.8, 1.0; 1.0, 1.1, 1.2; 1.5, 1.6, 1.7; 15, 16, 13.5];

if ms <= (4/alfa)
    fprintf('ms > 4/alpha not met')
    error
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
styleN = 'xo'

for j = 1:length(N)
    for k = 1:length(z)
        % Ao = sqrt(gammaBar*(2+z(k)^2))/(rc*z(k)*Hl);
        % for i = 1:length(gammaBar)
        %     [j k i]
        %     % Ganhos aleatórios
        %     Hf = gainAF(alfa(j),mu,ms,rc,Nc,-1e-3); % Alpha F
        %     Hp = PointError(z(k),Ao(i),Nc); % Pointing error
        %     % Ganho total
        %     Gain = (Hl(:).*Hf(:).*Hp(:)).';
            
        %     [flagOP,~] = find(Gain.^2 <= gamma_th);
        %     Pout(i) = sum(flagOP)/Nc;
        % end

        %%
        
        nj = N(j);
        [gammaBar_dB, Pb] = OP_analit(N(j), alfa, mu(1:nj), ms, bounds, points, z(k, 1:nj), gamma_th);
        [gammaBar_dB, P] = OP_asymptotic(N(j), alfa, mu(1:nj), ms, bounds, points, z(k, 1:nj), gamma_th);

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
        semilogy(gammaBar_dB, Pb, [colorz(k),styleN(nj)],...
                 gammaBar_dB, P,'k--',...
                 'linewidth',1.2)
        
        hold on
    end
end

axis([min(GBdB) max(GBdB) 1e-5 1])
tam_fonte = 11;
legend('FontSize', tam_fonte)
legend("$z_1="+num2str(z(1,1))+", z_2="+num2str(z(1,2))+"$", '', "$z_1="+num2str(z(2,1))+", z_2="+num2str(z(2,2))+"$", '', "$z_1="+num2str(z(3,1))+", z_2="+num2str(z(3,2))+"$", '', "Non-pointing errors",'Asymptotic', 'Location', 'southwest')
ax = gca;
ax.FontSize = 13;

%title('Outage Probability')
set(legend, 'Interpreter', 'latex')
ylabel("OP", 'FontSize', 14)
yticks([1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1e0])
xlabel("SNR (dB)", 'FontSize', 14)
grid on

%textbox com valores
dim = [0.15 0.40 0.2 0.2];
str = {"$\alpha ="+num2str(alfa)+"$", "$\mu_1 ="+num2str(mu(1))+"$", "$\mu_2 ="+num2str(mu(2))+"$" , "$m ="+num2str(ms)+"$", "$\gamma_{\rm th}="+num2str(gamma_th_dB)+"$ dB"};
annotation('textbox',dim,'interpreter','latex','String',str,'FitBoxToText','on', 'FontSize', tam_fonte);

%text arrow
X = [0.74,0.7];
Y = [0.84,0.71];
str = {"$N=\{" + num2str(N(1)) + "," + num2str(N(2)) + "\}$"};
annotation('textarrow', X, Y, 'String', str, 'interpreter', 'latex', 'FontSize', 13);
