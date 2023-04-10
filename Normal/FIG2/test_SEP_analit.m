clear all
close all
clc
 
%QPSK
theta = 2;
phi = 1;

alpha = [2, 6.5];  % non-linearity parameter
mu = [1, 1.7];     % number of multipath clusters
ms = [2.5, 150];    % shadowing parameter ms > 2 / alpha

L = 1;
U = db2pow(30);
bounds = [L U]; % 0 - 30 dB
N = 1e4;        % Number of points

hl = 1;
z = 2.3;

A_0 = 0.8;
colorstring = 'bgrm';

aux = 0;
cont = 0;
cont1 = 0;
cont2 = 0;
txt1 = [0,0,0,0];

figure()
for k = 1:length(alpha)
    if k == 1
        aux = 0;
    else
        aux = 2;
    end
    for j = 1:length(ms)
        for i = 1:length(mu)
            cont2 = cont2 + 1;
            
            [gammaBar_dB, Pb] = SEP_analit(theta, alpha(k), A_0, mu(i), ms(j), phi, bounds, N, hl, z);
            semilogy(gammaBar_dB, Pb, 'Color', colorstring(i+aux))
            %txt(i*2-1) = "BPSK";
            txt(cont2*2-1) = "\alpha = " + num2str(alpha(k)) + " z = " + num2str(z) + " \mu = " + num2str(mu(i));
            %cont = (i+aux)*2-1
            hold on
            
            [gammaBar_dB, P] = SEP_asymptotic(theta, alpha(k), A_0, mu(i), ms(j), phi, bounds, N, hl, z);
            semilogy(gammaBar_dB, P, 'Color', 'k', 'LineStyle', '--')
            txt(cont2*2) = "asymptote";
            %cont1 = (i+aux)*2
            legend(txt)
            title("QPSK: {m_s} = " + num2str(ms))
            
        end
    end
end

% %NBFSK
% theta = 2;
% phi = 1;
% colorstring = 'r';
% 
% [gammaBar_dB, Pb] = SEP_analit(theta, alpha, betaVar, mu, ms, phi, bounds, N);
% semilogy(gammaBar_dB, Pb, 'Color', colorstring)
% txt = "NBFSK";
% 
% [gammaBar_dB, P] = SEP_asymptotic(theta, alpha, betaVar, mu, ms, phi, bounds, N);
% semilogy(gammaBar_dB, P, 'Color', colorstring, 'LineStyle', '--')
% txt = [txt "asymptote"];
% legend(txt)

% alpha = 2;
% mu = 2;
% colorstring = 'g';
% 
% [gammaBar_dB, Pb] = SEP_analit(theta, alpha, betaVar, mu, ms, phi, bounds, N);
% semilogy(gammaBar_dB, Pb, 'Color', colorstring)
% txt = [txt ("\alpha = " + num2str(alpha) + " \mu = " + num2str(mu)) ];
% 
% [gammaBar_dB, P] = SEP_asymptotic(theta, alpha, betaVar, mu, ms, phi, bounds, N);
% semilogy(gammaBar_dB, P, 'Color', colorstring, 'LineStyle', '--')
% txt = [txt "asymptote"];
% legend(txt)

grid on
% xlabel()
ylim([1e-5, 1e0])
legend('Location','southwest')
