clear all
close all
clc

set(0,'defaulttextinterpreter','latex');

M = 2; % Ordem da constela��o M-QAM
N = 5e5; % N�mero de bits transmitidos

L = 0;  %db
U = 30; %db
points = 1e2;
bounds = [L U];

theta  = 1;
phi = 2;
SNR_dB = linspace(L,U,15);

%Parâmetros
mu = [1 , 1.7];
alfa = [6.5 , 2];
Hl = 1.0;
ms = [2.5 , 100];
z = 3;
Ao = 0.8;

ber = zeros(length(SNR_dB),length(M));
ser = zeros(length(SNR_dB),length(M));


colorstr1 = 'bg';
colorstr2 = 'cm';
v = 1;
for Alpha = 1:length(alfa)
    for k = 1:length(ms)
        for m = 1:length(mu)
            m
           for i=1:length(SNR_dB)
              Akmi = [Alpha k m i]
               [ser(i,m),ber(i,m)] = channel(M,SNR_dB(i),N,alfa(Alpha),mu(m),ms(k),z,Ao,Hl);
           end

            [gammaBar_dB, Pb] = SEP_analit(theta, alfa(Alpha), mu(m), ms(k), phi, bounds, points, z);
            [gammaBar_dB, P] = SEP_asymptotic(theta,alfa(Alpha), mu(m), ms(k), phi, bounds, points, z);

           
            if Alpha == 1
                figure(1)
                semilogy(SNR_dB,ser(:,m),'rx',...
                        gammaBar_dB, Pb, colorstr1(m),...
                        gammaBar_dB, P,'k--',...
                        'linewidth',1.5)
            end
            if Alpha == 2
                figure(1)
                semilogy(SNR_dB,ser(:,m),'rx',...
                        gammaBar_dB, Pb, colorstr2(m),...
                        gammaBar_dB, P,'k--',...
                        'linewidth',1.5)
            end
            hold on
        end
    end
end
grid on
%legend('Simulated', 'BPSK','','', 'QPSK', 'Asymptotic', 'Location', 'southwest')
legend('Simulated','alpha = 6.5  \mu = 0.3','Asymptotic','','alpha = 6.5  \mu = 2','','', 'alpha = 2  \mu = 0.3','','','alpha = 2  \mu = 2' , 'Location', 'southwest')

set(legend, 'Interpreter', 'latex')
%title("BPSK and QPSK")
ylim([1e-5, 1e0])
ylabel('SEP', 'FontSize', 14)
xlabel('SNR', 'FontSize', 14)

dim = [0.15 0.35 0.2 0.2];
str = {'$\alpha = 3$', '$\mu = 2$', '$m_s=5$', '$A_0=0.8$'};
annotation('textbox', dim, 'interpreter', 'latex', 'String', str, 'FitBoxToText', 'on');
