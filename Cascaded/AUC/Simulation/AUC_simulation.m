function [AUC] = AUC_simulation(Nc, N, params, gammaBar, u, N_auc, N_roc)
    % params is a matrix of params cascaded system
    % params = [Channel 1: alpha, mu, ms, z, rc, Hl;
    %           Channel 2: alpha, mu, ms, z, rc, Hl;
    %           ...
    %           Channel n: alpha, mu, ms, z, rc, Hl;]
    % cascaded_gain = gain_channel1 * gain_channel2 * ... * gain_channeln

    % N: number of channels
    % alpha: non-linearity
    % mu: number of multipaths
    % ms: shadowing
    % z: pointing error
    % Ao: pointing error param
    % rc: hat r (average power)
    % Nc: number of simulation points
    % Hl: path loss
    % u: time-bandwith
    % N_roc: points number in ROC curve (Receiver Operation Characteristics curve)
    % N_auc: points in AUC curve to calculate average AUC
    % gammaBar: SNR matriz per channel
    % gammaBar = [          , channel1, channel2, ..., channeln;
    %              gammaBar1,     y1.1,     y1.2, ...,     y1.n;
    %                    ...,      ...,      ..., ...,      ...;
    %             gammaBar15,    y15.1,    y15.2, ...,    y15.n;]

    Nu = u*2; % numero de amostras coletadas por radio com base em u (time-bandwith)
    
    x = randn(1, Nu);  % sinal transmitido
    Ptx = mean(x.^2); % potência do sinal transmitido
    x = repmat(x, N_roc, 1);

    energy_min = 0;
    energy_max = 300;

    energy = linspace(energy_min, energy_max, N_roc)';
    AUC = zeros(1, length(gammaBar));

    channels_gain = individual_gain(N, params, Nc, gammaBar);
    total_gain = cascaded_gain(channels_gain);
    Gain = total_gain(:,:,N);  % Gain: Nc(linhas) X gammaBar(colunas)

    % for i = 1:length(gammaBar)
    %     Gain(:, i) = Gain(:, i) ./ sqrt(mean(Gain(:, i).^2));
    % end
  

    AUC = zeros(1, length(gammaBar));
    for j = 1:length(gammaBar)
        j

        sum_AUC = 0;
        for k = 1:N_auc
            % k

            I = randi([0, 1], Nc, 1);   % variavel de bernoulli indicando se há ou não
                                        % usuário primário utilizando o canal
            % gerando ruído
            % sigma_ruido = sqrt((u*Ptx) ./ gammaBar(j, 2)); % desvio padrão ruído
            sigma_ruido = sqrt((u*Ptx) ./ gammaBar(j, N)); % desvio padrão ruído
            w = (sigma_ruido) .* randn(N_roc, Nu, Nc);

            cont_d = zeros(N_roc, Nc);
            cont_totald = 0;
            cont_fa = zeros(N_roc, Nc);
            cont_totalfa = 0;

            for i = 1:Nc
                % sinal recebido
                r = I(i).*Gain(i, j).*x + w(:,:,i);
                
                % energia do sinal recebido
                y = sum((r.^2),2)/(sigma_ruido^2);

                % contador de detecção e falso alarme baseado na energia
                cont_d(:, i)  = (abs(y) > energy) .* I(i);      % havia usuário e foi detectado
                cont_fa(:, i) = (abs(y) > energy) .* not(I(i)); % não havia usuário e foi detectado

                cont_totald = cont_totald + I(i);
                cont_totalfa = cont_totalfa + I(i);
            end

            % probabilidade de deteccao e falso alarme
            Pd  = ((sum(cont_d,2)) ./ cont_totald)';
            Pfa = ((sum(cont_fa,2)) ./ cont_totalfa)';

            % soma aucs para calcular a média
            sum_AUC = sum_AUC + abs(trapz([1 Pfa 0],[1 Pd 0]));
        end
        AUC(j) = sum_AUC/N_auc;
    end

end