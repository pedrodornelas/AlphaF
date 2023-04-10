clear all
close all
clc

set(0,'defaulttextinterpreter','latex');
p = [0,0,0; 0,0,0; 0,0,0];
q = [0,0,0; 0,0,0; 0,0,0];

M = [2 , 4]; % Ordem da constela��o M-QAM
N = 1e5; % N�mero de bits transmitidos

L = 0;  %db
U = 30; %db
points = 1e2;
bounds = [L U];

theta  = [1 , 2];
phi = [2 , 1];
SNR_dB = linspace(L,U,15);

%Parâmetros
mu = 2.1;
alfa = 3.5;
Hl = 1.0;
ms = 5;
z = [0.6, 1.1, 6.5];
Ao = 0.8;

ber = zeros(length(SNR_dB),length(M));
ser = zeros(length(SNR_dB),length(M));

colorstr = 'bg';
v = 1;
for k = 1:length(z)
    for m = 1:length(M)
        m
        for i=1:length(SNR_dB)
            kmi = [k m i]
            [ser(i,m),ber(i,m)] = channel(M(m),SNR_dB(i),N,alfa,mu,ms,z(k),Ao,Hl);
        end
        
        [gammaBar_dB, Pb] = SEP_analit(theta(m), alfa, mu, ms, phi(m), bounds, points, z(k));
        [gammaBar_dB, P] = SEP_asymptotic(theta(m),alfa, mu, ms, phi(m), bounds, points, z(k));
        
        figure(1)
        if m == 1
            p(: , k) = semilogy(SNR_dB,ser(:,m),'rx',...
                    gammaBar_dB, Pb, colorstr(m),...
                    gammaBar_dB, P,'k--',...
                    'linewidth',1.5)
        end
        
        if m == 2
            q(: , k) = semilogy(SNR_dB,ser(:,m),'rx',...
                    gammaBar_dB, Pb, colorstr(m),...
                    gammaBar_dB, P,'k--',...
                    'linewidth',1.5)
        end
        hold on
    end
end

grid on
legend([p(1) p(2) q(2) q(3)],'Simulated', 'BPSK','','', 'QPSK', 'Asymptotic', 'Location', 'southwest')
set(legend, 'Interpreter', 'latex')
%title("BPSK and QPSK")
ylim([1e-5, 1e0])
yticks([1e-5 1e-4 1e-3 1e-2 1e-1 1e0])
ylabel('SEP', 'FontSize', 14)
xlabel('SNR', 'FontSize', 14)

dim = [0.15 0.35 0.2 0.2];
str = {"$\alpha = "+num2str(alfa)+"$", "$\mu = "+num2str(mu)+"$", "$m_s= "+num2str(ms)+"$", "$A_0="+num2str(Ao)+"$"};
annotation('textbox', dim, 'interpreter', 'latex', 'String', str, 'FitBoxToText', 'on');
