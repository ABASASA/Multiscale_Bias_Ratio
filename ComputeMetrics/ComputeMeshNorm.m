% X, Y - are posistion vectors
% Assuming the space is a square
function [distance] = ComputeMeshNorm(X, Y, gridLimits)


sampleGrid = linspace(gridLimits(1), gridLimits(2), 1000);
n = length(sampleGrid);
XSpace = repmat(sampleGrid, 1, n);
XSpace = XSpace(:);
YSpace = repmat(sampleGrid, n, 1);
YSpace = YSpace(:);

distance = -inf;

for i = 1 : n^2
   
    tmpDist = min(sqrt((X - XSpace(i)).^2 + (Y - YSpace(i)).^2));
    distance = max(tmpDist, distance);
end


end