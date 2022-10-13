function efcp = fcpe(condition,vector)
% Função para a geração da FCP empírica (eFCP) das
% amostras contidas em "vector". Condition é um vetor
% linha contendo os valores que se deseja a eFCP.
% Para condition com um único elemento a função 
% não entra no loop.

efcp = zeros(1,length(condition)); % Alocação de Memória

% Laço para a montagem da eFCP
for i = 1:1:length(condition)
    % Vetor para realizar a comparação X <= x da eFCP
    vector_comp = ones(size(vector))*condition(i);
    % Valor da eFCP para X <= x (x=condition(i))
    efcp(i) = length(find(vector<=vector_comp))/length(vector);
end

end
