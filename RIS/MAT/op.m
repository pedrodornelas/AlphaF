clear all
close all

load('OP_L_1_gamma_th_5_dB.mat');
OP1=OP;
OP_Asymp1=OP_Asymp;

load('OP_L_2_gamma_th_5_dB.mat');
OP2=OP;
OP_Asymp2=OP_Asymp;

load('OP_L_3_gamma_th_5_dB.mat');
OP3=OP;
OP_Asymp3=OP_Asymp;

load('OP_L_4_gamma_th_5_dB.mat');
OP4=OP;
OP_Asymp4=OP_Asymp;

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
