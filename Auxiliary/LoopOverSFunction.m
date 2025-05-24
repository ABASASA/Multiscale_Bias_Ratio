function [ss] = LoopOverSFunction(s, X, nDimOut)

ss = zeros(size(X,1), nDimOut);

for i = 1 : length(X)
    ss(i, :) = s(X(i, :));
end

end