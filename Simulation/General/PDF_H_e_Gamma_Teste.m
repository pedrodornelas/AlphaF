clear all
clc

% Piece of code to validate the tehoretical expression
% Brito - 28/09/2022

% Parâmetros da Distribuição Alfa F
alfa = 3.5;
mu = 2;
ms = 5;
rc = 1;
Nc = 1e6;

% Parametros da distribuição do erro de apontamento
z = 0.80;
Ao = 1.0;

% Perda de percurso
Hl = 1.00;

% Ganhos aleatórios
Hf = gainAF(alfa,mu,ms,rc,Nc,1e-3); % Alpha F
Hp = PointError(z,Ao,Nc); % Pointing error

% SNR média
gb = ((rc*Ao*z*Hl)^2)/(2+z^2);

% Ganho total
Gain = (Hl(:).*Hf(:).*Hp(:)).';

% Parâmetros da PDF teórica
% % Parâmetros de entrada da MeijerG
m = 2; n = 1;
p = 2; q = 2;
a = [1-ms 1+z*z/alfa];
b = [mu z*z/alfa];

%%

Psi = mu/(ms-1);
% PDF Gain
pdfH =@(h) z*z/gamma(mu)/gamma(ms)./h.*...
          meijerG(a(1:n),a(n+1:end),b(1:m),b(m+1:end),Psi*(h/rc/Hl/Ao).^alfa);

% PDF SNR
pdfG =@(g) z*z/gamma(mu)/gamma(ms)/2./g.*...
          meijerG(a(1:n),a(n+1:end),...
                  b(1:m),b(m+1:end),...
                  Psi*(z*sqrt(g/gb/(2+z^2))).^alfa);


% % Parâmetros de entrada da MeijerG
m = 2; n = 2;
p = 3; q = 3;
a = [1-ms 1 1+z*z/alfa];
b = [mu z*z/alfa 0];
% PDF SNR
cdfG =@(g) z*z/alfa/gamma(mu)/gamma(ms).*...
          meijerG(a(1:n),a(n+1:end),...
                  b(1:m),b(m+1:end),...
                  Psi*(z*sqrt(g/gb/(2+z^2))).^alfa);

% Histograma normalizado -- Gain
[fx,x] = histnorm(Gain,1e2);
% Histograma normalizado -- SNR
[fy,y] = histnorm(Gain.^2,1.5e2);
%%

[G,CG] = ecdf(Gain.^2);

% Gráficos
figure(1)
plot(x,pdfH(x),...
     x,fx,'--',...
     'linewidth',1.5)

figure(2)
plot(y,pdfG(y),...
     y,fy,'--',...
     'linewidth',1.5)
%%
clc

GG = linspace(1e-5,4,100);
CG = fcpe(GG,Gain.^2);

figure(3)
plot(GG,cdfG(GG),...
     GG,CG,'--',...
     'linewidth',1.5)
axis([0 4 0 1])

%%
% Segundo Momento amostral de Hf (ou SNR média amostral)
Rc = mean(Hf.^2)
% Rc média teórica
rc^2
% Segundo Momento amostral de H (ou SNR média amostral)
gba = mean(Gain.^2)

% SNR média teórica
gb = ((rc*Ao*z*Hl)^2)/(2+z^2)


%%

figure(123564231)
ecdf(Gain.^2)