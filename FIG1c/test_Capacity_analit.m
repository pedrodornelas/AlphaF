clear all
close all
clc

alpha = 2.7;  % non-linearity parameter
mu = 3;     % number of multipath clusters
ms = [3.1, 10.5];    % shadowing parameter ms > 2 / alpha

L = 1;
U = db2pow(30);
bounds = [L U]; % 0 - 30 dB
N = 1e4;        % Number of points

hl = 1;
z = [0.8, 1.5, 6.7];

gamma_th = 5;

colorstring = 'bgrm';
A_0 = 0.8;

figure()
for j = 1:length(ms)
    for i = 1:length(z)
        
        [gammaBar_dB, Pb] = Capacity_analit(alpha, A_0, mu, ms(j), bounds, N, hl, z(i));
        plot(gammaBar_dB, Pb, 'Color', colorstring(i))
        %txt(i*2-1) = "BPSK";
        txt(i*2-1) = "\alpha = " + num2str(alpha) + " \mu = " + num2str(mu) + " z = " + num2str(z(i));
        hold on
        
        [gammaBar_dB, P] = Capacity_asymptotic(alpha, A_0, mu, ms(j), bounds, N, hl, z(i));
        plot(gammaBar_dB, P, 'Color', colorstring(i), 'LineStyle', '--')
        txt(i*2) = "asymptote";
        legend(txt)
        title("ms = " + num2str(ms))
        
    end
end

% %NBFSK
% theta = 2;
% phi = 1;
% colorstring = 'r';
% 
% [gammaBar_dB, Pb] = Capacity_analit(theta, alpha, betaVar, mu, ms, phi, bounds, N);
% semilogy(gammaBar_dB, Pb, 'Color', colorstring)
% txt = "NBFSK";
% 
% [gammaBar_dB, P] = Capacity_asymptotic(theta, alpha, betaVar, mu, ms, phi, bounds, N);
% semilogy(gammaBar_dB, P, 'Color', colorstring, 'LineStyle', '--')
% txt = [txt "asymptote"];
% legend(txt)

% alpha = 2;
% mu = 2;
% colorstring = 'g';
% 
% [gammaBar_dB, Pb] = Capacity_analit(theta, alpha, betaVar, mu, ms, phi, bounds, N);
% semilogy(gammaBar_dB, Pb, 'Color', colorstring)
% txt = [txt ("\alpha = " + num2str(alpha) + " \mu = " + num2str(mu)) ];
% 
% [gammaBar_dB, P] = Capacity_asymptotic(theta, alpha, betaVar, mu, ms, phi, bounds, N);
% semilogy(gammaBar_dB, P, 'Color', colorstring, 'LineStyle', '--')
% txt = [txt "asymptote"];
% legend(txt)

grid on
% xlabel()
ylim([0, 10])
legend('Location', 'northwest')
