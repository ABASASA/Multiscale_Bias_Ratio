function [rates] = CreateRate(mu, nLevels, nDim)

rates = ones(nLevels, 1) ./ mu.^([(nLevels-1):-1:0]' ./ nDim);

end