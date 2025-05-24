% X, Y - are posistion vectors
% Assuming the space is a square
function [distance] = ComputeMeshNormRn(X, gridLimits, nDim)

if( nDim ~= size(X,2))
    error('Code Error: Dimension probelm');
end

sampleGrid = linspace(gridLimits(1), gridLimits(2), 100);

nPointsInAxis = length(sampleGrid);

% First stage is to  built the list of point to check iteratavlly
points = zeros(nPointsInAxis.^nDim, nDim);

if nDim == 1
    points = sampleGrid';
else
    r = ones(1, nDim);
    for i = 1 : nDim
       r(i) =  nPointsInAxis;
       tmpp = repmat(sampleGrid, r);
       points(:, i) = tmpp(:);
       r(i) = 1;
    end
end

distance = -inf;

for i = 1 : nPointsInAxis^nDim
    tmpPoint = points(i,:);
    summ = 0;
    for j = 1 : nDim
        summ = summ + (X(:, j) - tmpPoint(j)).^2;
    end
    tmpDist = min(sqrt(summ));
    distance = max(tmpDist, distance);
end


end