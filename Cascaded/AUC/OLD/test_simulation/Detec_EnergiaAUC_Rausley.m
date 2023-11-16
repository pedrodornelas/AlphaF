clc
clear all
close all

%% Simulado
N_runs = 10^4; % n�mero de eventos de Monte Carlo para o c�lculo de cada valor de AUC
N_limiar = 50; % n�meros de pontos na ROC
N = 20; % n�mero de amostras coletadas por radio
u=N/2;
N_auc = 100;%n�mero de valores de AUCs que ser� usado para o c�lculo da AUC m�dia

SNR = 8; %em dB
SNR_mean = 10.^(SNR./10);

x = randn(1,N); %sinal transmitido determin�stico 
Ptx = mean(x.^2);
sigma_ruido = sqrt((N.*Ptx)./(2.*SNR_mean));
x = repmat(x,N_limiar,1);

l = N_runs;

limiar_min = 0;
limiar_max = 300;

limiar = (linspace(limiar_min,limiar_max,N_limiar))';

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
    
               I = randi([0 1]);
               w = (sigma_ruido).*randn(N_limiar,N);
    
               r = I.*h.*x + w;
    
               y(:,i) = (sum((r.^2),2)/(sigma_ruido^2));
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

AUC = sum(AUC_S)./N_auc
