function [errorSource, pointsTest] = ComputeErrorSPD(fEstimate, func,...
                        gridLimits, nDim, flag, outliersPos, testDist, flagAllLevels)

if ~exist('flag', 'var')
    flag = false;
end

if ~exist('outliersPos','var')
    outliersPos = [];
end
  
if ~exist('testDist','var')
    testDist = 0.01;
end

if ~exist('flagAllLevels','var')
    flagAllLevels = true;
end


if flag
    figure;
end

[pointsTest, ~,~] = SampleRnfunctionEqui(func, gridLimits, testDist, nDim);
nDimOut = size(pointsTest, 2) - nDim;

% errorSource = zeros(size(pointsTest,1), length(fEstimate)); 
errorScales = zeros(size(pointsTest,1), length(fEstimate));

tmp = zeros(size(pointsTest,1), nDimOut);
sizeSide = sqrt(size(pointsTest,1));
counter = 1;
tmpOld = [];

if flagAllLevels
    startLevel = 1;
    errorSource = zeros(size(pointsTest,1), length(fEstimate)); 
else
    startLevel = length(fEstimate);
    errorSource = zeros(size(pointsTest,1), 1); 
end

for iLevel = startLevel : length(fEstimate)
    tmpFucntion = fEstimate{iLevel};
    
    parfor iPoint = 1 : size(pointsTest,1)
        tmp(iPoint, :) = tmpFucntion(pointsTest(iPoint,1:nDim));
    end
    if flagAllLevels

        errorSource(:, iLevel) = GeodesicDistanceSPD(pointsTest(:,nDim + 1 : end), tmp);
    else
        errorSource = GeodesicDistanceSPD(pointsTest(:,nDim + 1 : end), tmp);
    end
    
    if (iLevel ~= 1 && flagAllLevels)
        errorScales(:, iLevel - 1) = GeodesicDistanceSPD(oldTmp, tmp);
    end
    oldTmp = tmp;
    

end

end

