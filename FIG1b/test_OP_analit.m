clear all
close all
clc

alpha = [2.0, 3.5];  % non-linearity parameter
mu = 3;     % number of multipath clusters
ms = 1.3;    % shadowing parameter ms > 2 / alpha

L = 1;
U = db2pow(30);
bounds = [L U]; % 0 - 30 dB
N = 1e4;        % Number of points

hl = 1;
z = [0.8, 1.5, 6.7]; % near, moderate, stroug pointing error

gamma_th_dB = 5;
gamma_th = db2pow(gamma_th_dB);

A_0 = 0.8;
colorstring = 'bgrm';

lstyle = '';

figure()
for j = 1:length(alpha)
    %j=1 -> alpha = 2.0
    %j=2 -> alpha = 3.5
    if j == 1
        lstyle = '-';
    else
        lstyle = '-';
    end
    for i = 1:length(z)
        
        [gammaBar_dB, Pb] = OP_analit(alpha(j), A_0, mu, ms, bounds, N, hl, z(i), gamma_th);
        semilogy(gammaBar_dB, Pb, 'Color', colorstring(i), 'LineStyle', lstyle)
        %txt(i*2-1) = "BPSK";
        txt(i*2-1) = "\alpha = " + num2str(alpha(j)) + " \mu = " + num2str(mu) + " z = " + num2str(z(i));
        hold on
        
        [gammaBar_dB, P] = OP_asymptotic(alpha(j), A_0, mu, ms, bounds, N, hl, z(i), gamma_th);
        semilogy(gammaBar_dB, P, 'Color', colorstring(i), 'LineStyle', '--')
        txt(i*2) = "asymptote";
        legend(txt)
        title("\alpha = " + num2str(alpha))
        
    end
end

% %NBFSK
% theta = 2;
% phi = 1;
% colorstring = 'r';
% 
% [gammaBar_dB, Pb] = OP_analit(theta, alpha, betaVar, mu, ms, phi, bounds, N);
% semilogy(gammaBar_dB, Pb, 'Color', colorstring)
% txt = "NBFSK";
% 
% [gammaBar_dB, P] = OP_asymptotic(theta, alpha, betaVar, mu, ms, phi, bounds, N);
% semilogy(gammaBar_dB, P, 'Color', colorstring, 'LineStyle', '--')
% txt = [txt "asymptote"];
% legend(txt)

% alpha = 2;
% mu = 2;
% colorstring = 'g';
% 
% [gammaBar_dB, Pb] = OP_analit(theta, alpha, betaVar, mu, ms, phi, bounds, N);
% semilogy(gammaBar_dB, Pb, 'Color', colorstring)
% txt = [txt ("\alpha = " + num2str(alpha) + " \mu = " + num2str(mu)) ];
% 
% [gammaBar_dB, P] = OP_asymptotic(theta, alpha, betaVar, mu, ms, phi, bounds, N);
% semilogy(gammaBar_dB, P, 'Color', colorstring, 'LineStyle', '--')
% txt = [txt "asymptote"];
% legend(txt)

grid on
% xlabel()
ylim([1e-5, 1e0])
legend('Location','southwest')
