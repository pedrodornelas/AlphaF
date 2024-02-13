clear all;
% close all;

set(0,'defaulttextinterpreter','latex');
set(groot,'defaultAxesTickLabelInterpreter','latex');

load('BEP.mat');

figure()
colorz = 'brgmp';
cont = 1;
h = [];

for i = 1:length(L)
    for j = 1:length(z)
        h(cont)   = semilogy(gamma_bar_dB, BEP(:, i, j), colorz(i), 'linewidth', 1.2); hold on;
        h(cont+1) = semilogy(gamma_bar_dB, BEP_asy(:, i, j), 'k--', 'linewidth', 1); hold on;
        cont = cont+2;
    end
end


axis([min(gamma_bar_dB) max(gamma_bar_dB) 1e-5 1])
tam_fonte = 11;
legend('FontSize', tam_fonte, 'Location', 'southeast')
legend([h(1), h(5), h(9), h(13), h(14)], {"$L="+num2str(L(1))+"$", "$L="+num2str(L(2))+"$", "$L="+num2str(L(3))+"$", "$L="+num2str(L(4))+"$", "Asymptotic"});
% legend([h(1), h(5), h(9), h(10)], {"$L="+num2str(L(1))+"$", "$L="+num2str(L(2))+"$", "$L="+num2str(L(3))+"$", "Asymptotic"});
ax = gca;
ax.FontSize = 13;

%title('Outage Probability')
set(legend, 'Interpreter', 'latex')
ylabel("BEP", 'FontSize', 14)
yticks([1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1e0])
xlabel("SNR (dB)", 'FontSize', 14)
grid on

%textbox com valores
% dim = [0.148571428571428,0.333971669950265,0.242555404298264,0.226504520525927];
% str = {"$\alpha_1 ="+num2str(alpha(1))+", \alpha_2 ="+num2str(alpha(2))+"$",...
%        "$\mu_1 ="+num2str(mu(1))+", \mu_2 ="+num2str(mu(2))+"$",...
%        "$m_1 ="+num2str(ms(1))+", m_2 ="+num2str(ms(2))+"$",...
%        "$\rho="+num2str(p)+"$ (BPSK)"};
% annotation('textbox',dim,'interpreter','latex','String',str,'FitBoxToText','on', 'FontSize', tam_fonte);