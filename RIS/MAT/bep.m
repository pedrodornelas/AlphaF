clear all
close all

load('BEP_L_1_p_1.mat');
OP1=BEP;
OP_Asymp1=BEP_Asymp;

load('BEP_L_2_p_1.mat');
OP2=BEP;
OP_Asymp2=BEP_Asymp;

load('BEP_L_3_p_1.mat');
OP3=BEP;
OP_Asymp3=BEP_Asymp;

load('BEP_L_4_p_1.mat');
OP4=BEP;
OP_Asymp4=BEP_Asymp;

figure()
colorstring = 'bgrm';
semilogy(gamma_bar_dB, OP1,'Color','b','LineStyle','-')
hold on
semilogy(gamma_bar_dB, OP_Asymp1,'Color','k','LineStyle','--')

hold on
semilogy(gamma_bar_dB, OP2,'Color','g','LineStyle','-')
hold on
semilogy(gamma_bar_dB, OP_Asymp2,'Color','k','LineStyle','--')

hold on
semilogy(gamma_bar_dB, OP3,'Color','r','LineStyle','-')
hold on
semilogy(gamma_bar_dB, OP_Asymp3,'Color','k','LineStyle','--')

hold on
semilogy(gamma_bar_dB, OP4,'Color','m','LineStyle','-')
hold on
semilogy(gamma_bar_dB, OP_Asymp4,'Color','k','LineStyle','--')


grid on
