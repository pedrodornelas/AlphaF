function point = PointError(z,Ao,NF)

PDF =@(y) (z*z*y.^(z*z-1))/Ao^(z*z);


y = PDF(0:0.01:Ao);
mx = real(max(y));
if mx > 100
    mx = 100;
end

gain = zeros(NF,1);
g = [];
while length(g) <= NF
    x1 = random('unif',0,Ao,[1 NF]);
    y = random('unif',0,mx,[1 NF]);
    [~,x] = find(y <= PDF(x1));
        g = cat(2,g,x1(x));
end


point(:,1) = g(1:NF); % Limitando o tamanho do vetor 

% [fx,x] = histnorm(point(:),1e2);
% e = 0:0.001:max(point(:));
% 
% figure(132)
% plot(e,PDF(e),'r',...
%      x,fx,'b:',...
%      'linewidth',1.5)

end
