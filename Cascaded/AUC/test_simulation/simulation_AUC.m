function [AUC] = simulation_AUC(Nc, gammaBar)

%% Simulado
N_runs = Nc; % numero de eventos de Monte Carlo para o calculo de cada valor de AUC
N_limiar = 1000; % numeros de pontos na ROC
N = 8; % numero de amostras coletadas por radio
u = N/2;
N_auc = 100;%numero de valores de AUCs que sera usado para o calculo da AUC media

% L = -10;  %dB
% U = 30;   %dB
% gammaBar_dB = linspace(L, U, 41);
% gammaBar = 10.^(gammaBar_dB./10);

% SNR = 8; %em dB
% SNR_mean = 10.^(SNR./10);

x = randn(1,N);   %sinal transmitido deterministico 
Ptx = mean(x.^2);
x = repmat(x,N_limiar,1);

l = N_runs;

limiar_min = 0;
limiar_max = 300;

limiar = (linspace(limiar_min,limiar_max,N_limiar))';
AUC = zeros(1, length(gammaBar));

for j = 1:length(gammaBar)
     j

     sigma_ruido = sqrt((N.*Ptx) ./ (2.*gammaBar(j)));
     AUC_S = zeros(N_auc,1);

     %for v=1:N_pontos     
     for f=1:N_auc
          w = zeros(N_limiar,N);
          r = zeros(N_limiar,N);

          cont_d = zeros(N_limiar,N_runs);
          cont_fa = zeros(N_limiar,N_runs);
          cont_totald = 0;
          cont_totalfa = 0;

          for i=1:N_runs
               % h = 1; %canal AWGN
               h = sqrt(1/2.*randn.^2+1/2.*randn.^2); %canal Rayleigh
               % h = h./mean(Pot_h^2);

               I = randi([0 1]);
               w = (sigma_ruido).*randn(N_limiar,N);
     
               % r = I.*h.*x + w;
               r = I.*h.*x + w;
               % size(r)
               
               y(:,i) = (sum((r.^2),2)/(sigma_ruido^2));
               % y(:,i) = (sum((r.^2),2));
               cont_d(:,i)  = (abs(y(:,i)) > limiar) .* I;
               cont_fa(:,i) = (abs(y(:,i)) > limiar) .* not(I);
          
               cont_totald = cont_totald + I;
               cont_totalfa = cont_totalfa + not(I);
          end
          
          Pd = ((sum(cont_d,2))./cont_totald)';
          Pfa = ((sum(cont_fa,2))./cont_totalfa)';    

          AUC_S(f)= abs(trapz([1 Pfa 0],[1 Pd 0]));
     end
     %end

     AUC(j) = sum(AUC_S)./N_auc;
end

% figure(1)
% plot(gammaBar_dB, AUC, 'rx')
end
