function efcp = fcpe(condition,vector)
% Fun��o para a gera��o da FCP emp�rica (eFCP) das
% amostras contidas em "vector". Condition � um vetor
% linha contendo os valores que se deseja a eFCP.
% Para condition com um �nico elemento a fun��o 
% n�o entra no loop.

efcp = zeros(1,length(condition)); % Aloca��o de Mem�ria

% La�o para a montagem da eFCP
for i = 1:1:length(condition)
    % Vetor para realizar a compara��o X <= x da eFCP
    vector_comp = ones(size(vector))*condition(i);
    % Valor da eFCP para X <= x (x=condition(i))
    efcp(i) = length(find(vector<=vector_comp))/length(vector);
end

end
