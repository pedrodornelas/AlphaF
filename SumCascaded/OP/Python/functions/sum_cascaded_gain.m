% function Gain = sum_cascaded_gain(L, N, gain_cascaded)
function Gain = sum_cascaded_gain(L, N, Nc, points, gammaBar, params)

    % Nc = size(gain_cascaded_L, 1);
    % gammaBar_points = size(gain_cascaded_L, 2);
    % cascate = size(gain_cascaded_L, 3);

    % Gain = zeros(Nc, gammaBar_points);
    % for l = 1:L
    %     for i = 1:gammaBar_points
    %         % Gain(:, i) = Gain(:, i) + gain_cascaded(:, i, N, l);
    %         Gain(:, i, l) = l .* gain_cascaded(:, i, N);
    %     end
    % end

    Gain = zeros(Nc, points, length(L));
    % ganho em cada cascata das L cascatas existentes
    gain_cascaded_L = zeros(Nc, points, length(L));

    for i = 1:length(L)
        L(i)
        gain_channels = individual_gain(max(N), params, Nc, gammaBar);
        gain_cascaded = cascaded_gain(gain_channels);
        gain_cascaded_L(:, :, i) = gain_cascaded(:, :, N);
    end

    for i = 1:length(L)
        if i == 1
            Gain(:,:,i) = gain_cascaded_L(:,:,i);
        else
            Gain(:,:,i) = Gain(:,:,(i-1)) + gain_cascaded_L(:,:,i);
        end
    end
    
end