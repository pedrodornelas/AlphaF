clear all
close all
clc
 
%BPSK
theta = 1;
phi = 2;

alpha = 3.5;  % non-linearity parameter
mu = 3;       % number of multipath clusters
ms = 5;       % shadowing parameter ms > 2 / alpha

L = 1;
U = db2pow(30);
bounds = [L U]; % 0 - 30 dB
N = 1e4;        % Number of points

z = [0.8, 1.5, 6.7]; % stroug, moderate and near pointing error

gamma_th_dB = 5; %in dB
gamma_th = db2pow(gamma_th_dB);

colorstring = 'bgrm';

figure()
for i = 1:length(z)
    
    [gammaBar_dB, Pb] = OP_analit(alpha, mu, ms, bounds, N, z(i), gamma_th);
    semilogy(gammaBar_dB, Pb, 'Color', colorstring(i))
    %txt(i*2-1) = "BPSK";
    txt(i*2-1) = "\alpha = " + num2str(alpha) + " \mu = " + num2str(mu) + " z = " + num2str(z(i));
    i*2-1
    hold on
     
    [gammaBar_dB, P] = OP_asymptotic(alpha, mu, ms, bounds, N, z(i), gamma_th);
    semilogy(gammaBar_dB, P, 'Color', colorstring(i), 'LineStyle', '--')
    txt(i*2) = "asymptote";
    i*2
    legend(txt)
    title("z = " + num2str(z))
    
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
ylim([1e-8, 1e0])
legend('Location', 'southwest')