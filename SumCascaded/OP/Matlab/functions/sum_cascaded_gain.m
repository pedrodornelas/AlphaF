function Gain = sum_cascaded_gain(L, N, gain_cascaded)

    Nc = size(gain_cascaded, 1);
    gammaBar_points = size(gain_cascaded, 2);
    cascate = size(gain_cascaded, 3);

    Gain = ones(Nc, gammaBar_points, L);
    for l = 1:L
        for i = 1:gammaBar_points
            Gain(:, i, l) = l .* gain_cascaded(:, i, N);
        end
    end
    
end