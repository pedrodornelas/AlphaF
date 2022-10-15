

clc
clear all

M = [2 , 4]; % Ordem da constela��o M-QAM
N = 5e5; % N�mero de bits transmitidos

L = -30;  %db
U = 30; %db
points = 1e2;
bounds = [L U];

theta  = [1 , 2];
phi = [2 , 1];
SNR_dB = linspace(0,30,15);
mu = 2;
alfa = 3.0;
Hl = 1.0;
ms = 5;
z = 0.8;
Ao = 0.8;

ber = zeros(length(SNR_dB),length(M));
ser = zeros(length(SNR_dB),length(M));
v = 1;
for m = 1:length(M)
    m;
    for i=1:length(SNR_dB)
        mi = [m i]
        [ser(i,m),ber(i,m)] = channel(M(m),SNR_dB(i),N,alfa,mu,ms,z,Ao,Hl);
    end
    [gammaBar_dB, Pb] = SEP_analit(theta(m), alfa, mu, ms, phi(m), bounds, points, z);
    [gammaBar_dB, P] = SEP_asymptotic(theta(m),alfa, mu, ms, phi(m), bounds, points, z);
    figure(1)
    semilogy(SNR_dB,ber,'rx',...
             gammaBar_dB, P,'k--',...
             gammaBar_dB, Pb,'b',...
             'linewidth',1.5)
    hold on
end

grid on
legend('Simulated', '','Asymptotic', 'Theoretical', 'Location', 'southwest')
title("BPSK and QPSK")
ylim([1e-5, 1e0])