function [g] = gainAF(x,mu,ms,rc,N,tol)
% Função para geração de variáveis aleatórias I.I.D. segundo uma
% distribuição alpha-F de parâmetros:
% alpha = x
% mu = mu
% ms = ms
% rc = rc
% Com número de amostras N e com suporte contendo 100(1-tol)% da área
% abaixo da curva da FDP -- Pelo bem da tua alma, nunca coloques tol = 0.
% Valores razoáveis de tol estão acima de 5e-3.

% Casos particulares:
% x = 2 -- Fisher-Snedecor
% ms -> infinito -- Alpha-Mu (valore muito grandes de ms podem gerar
% resultados não confiáveis devido à limitacões no rastreio do suporte,
% além de demora na geração dos dados).

% Chute inicial para o suporte
mn = 0; % inicio do intervalo
mx = 5; % fim do intervalo

% Fading - Envelope gain
O = (ms-1)/mu; % Eq. 1
C = x*(O)^ms*rc^(x*ms)/beta(mu,ms); % Eq. 2
f =@(g) C*g.^(x*mu-1).*(g.^x+O*rc^x).^(-mu-ms); % Eq. 3


% Trecho utilizado para tentar rastrear o suporte da função e melhorar a
% aderência das amostras à distribuição teórica (3)

% Vetor vec utilizado para mera vizualização das curvas teóricas e
% simuladas -- Remover comentário do trecho final do código
vec = linspace(mn,mx,1e5);
% Gargalo do código T_T. Preciso desenvolver um método melhor
while trapz(vec,f(vec)) < abs(1-tol) 
    vec = linspace(mn,mx,1e3);
    mx = 5 + mx;
end
vec = linspace(mn,mx,1e3);

% Definição do limiar superior do método de aceitação-rejeição -- Limiar
% fixo. Caso máximo superior a 10, trucar valor em pi.
mg = real(max(f(vec)));
if mg > 10
    mg = pi;
end


% Método da aceitação-rejeição
g = [];
while length(g) < N
    x1 = random('unif',mn,mx,[1 N]);
    y1 = random('unif',0,mg,[1 N]);
    [~,xg] = find(y1 <= f(x1));
        g = cat(2,g,x1(xg));
end

% Seleção de um vetor de comprimento N
g = g(1:N);


% Visualização da aderência dos dados -- remover comentaro das partes abaixo
% [fx,x] = histnorm(g,1.5e2);
% figure(2)
% plot(x,fx,'rx',...
%      vec,f(vec),'b',...
%      'linewidth',1.5)
end