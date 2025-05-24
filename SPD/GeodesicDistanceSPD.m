function [dist] = GeodesicDistanceSPD(X,Y)
if size(X, 1) == 1 | size(X, 2) == 1
    [dist] = GeodesicDistanceSPDOne(X, Y);
else 
    [dist] = GeodesicDistanceSPDMultiple(X, Y);
end
end

function [dist] = GeodesicDistanceSPDMultiple(Xs, Ys)
numMat = size(Xs, 1);
numEl = size(Xs, 2);
sqNumel = sqrt(numEl);

dist = zeros(numMat, 1);

for ii = 1 : numMat
    YMat = reshape(Ys(ii,:), [sqNumel, sqNumel]);
    XMat = reshape(Xs(ii,:), [sqNumel, sqNumel]);

    Xsqinv = sqrtm(inv(XMat));

    dist(ii) = norm(logm(Xsqinv * YMat * Xsqinv), 'fro');
end
end


function [dist] = GeodesicDistanceSPDOne(X, Y)

numEl = length(X);
sqNumel = sqrt(numEl);

YMat = reshape(Y, [sqNumel, sqNumel]);
XMat = reshape(X, [sqNumel, sqNumel]);

Xsqinv = sqrtm(inv(XMat));

dist = norm(logm(Xsqinv * YMat * Xsqinv), 'fro');

end