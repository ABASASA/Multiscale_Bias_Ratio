% Assume that the data is passed as vectors
function [results] = ExpSPD(Xs, Y)

nMats = size(Xs, 1);
numEl = size(Xs, 2);
sqNumel = sqrt(numEl);

results = zeros(nMats, numEl);

YMat = reshape(Y, [sqNumel, sqNumel]);
Ysq = sqrtm(YMat);
Ysqinv = inv(Ysq);
% Computing the Log for every matrxi in Xs.
for ii = 1 : nMats
    X = reshape(Xs(ii,: ), [sqNumel, sqNumel]);
    XexpY = Ysq * expm(Ysqinv * X * Ysqinv) * Ysq;
    results(ii, :) = reshape(XexpY, [numEl, 1]);
end


end

