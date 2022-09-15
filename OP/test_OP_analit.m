clear all
close all
clc
 
%BPSK
theta = 1;
phi = 2;

alpha = 3.5;  % non-linearity parameter
mu = 3;     % number of multipath clusters
ms = 5;    % shadowing parameter ms > 2 / alpha

L = 1;
U = db2pow(30);
bounds = [L U]; % 0 - 30 dB
N = 1e4;        % Number of points

hl = 1;
z = 20;

gamma_th_dB = 5;
gamma_th = db2pow(gamma_th_dB);

%betaVar = [1/3, 1, 2];      % AWGN
betaVar = [0.8, 1, 2];
colorstring = 'bgrm';

figure()
for i = 1:length(betaVar)
    
    [gammaBar_dB, Pb] = OP_analit(alpha, betaVar(i), mu, ms, bounds, N, hl, z, gamma_th);
    semilogy(gammaBar_dB, Pb, 'Color', colorstring(i))
    %txt(i*2-1) = "BPSK";
    txt(i*2-1) = "\alpha = " + num2str(alpha) + " \mu = " + num2str(mu) + "\beta = " + num2str(betaVar(i));
    i*2-1
    hold on
     
    [gammaBar_dB, P] = OP_asymptotic(alpha, betaVar(i), mu, ms, bounds, N, hl, z, gamma_th);
    semilogy(gammaBar_dB, P, 'Color', colorstring(i), 'LineStyle', '--')
    txt(i*2) = "asymptote";
    i*2
    legend(txt)
    title("\beta = " + num2str(betaVar))
    
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
