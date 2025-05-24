% Assume that the data is passed as vectors
function [results] = LogSPD(Xs, Y)
if size(Xs, 1) == size(Y, 1)
    [results] = LogSPDMultipleBases(Xs, Y);
else
    [results] = LogSPDOneBase(Xs, Y);
end
end

function [results] = LogSPDOneBase(Xs, Y)
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
    A = Ysqinv * X * Ysqinv;
%     if min(cond(A)) <= 0
%         1;
%     end
    XlogY = Ysq * logm(A) * Ysq;
    results(ii, :) = reshape(XlogY, [numEl, 1]);
end

end

function [results] = LogSPDMultipleBases(Xs, Ys)
nMats = size(Xs, 1);
numEl = size(Xs, 2);
sqNumel = sqrt(numEl);

results = zeros(nMats, numEl);
for ii = 1 : nMats
    YMat = reshape(Ys(ii,:), [sqNumel, sqNumel]);
    Ysq = sqrtm(YMat);
    Ysqinv = inv(Ysq);
    X = reshape(Xs(ii,: ), [sqNumel, sqNumel]);
    XlogY = Ysq * logm(Ysqinv * X * Ysqinv) * Ysq;
    results(ii, :) = reshape(XlogY, [numEl, 1]);
end  
    
end