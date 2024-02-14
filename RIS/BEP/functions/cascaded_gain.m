function Gain = cascaded_gain(gain_channels)

Nc = size(gain_channels, 1);
gammaBar_points = size(gain_channels, 2);
channels = size(gain_channels, 3);

Gain = ones(Nc, gammaBar_points, channels);
for c = 1:channels
    for i = 1:gammaBar_points
        if c == 1
            Gain(:, i, c) = gain_channels(:, i, c);
        else
            % Gain(:, i, c) = Gain(:, i, (c-1)) .* gain_channels(:, i, c);
            Gain(:, i, c) = Gain(:, i, (c-1)) .* gain_channels(:, i, c);
        end
    end
end

end