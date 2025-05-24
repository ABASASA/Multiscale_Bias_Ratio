function [errorMean, errorMedian, errorMax ] =...
        ComputeErrorAccording3MeasuremntsByLevel(bankOfEstimations, fucnTest,...
                                          gridLimitsTest, distTest, nRep, nDim, spaceData)


nLevels = length(bankOfEstimations{1});
errorMean = zeros(nRep, nLevels);
errorMedian = zeros(nRep, nLevels);
errorMax = zeros(nRep, nLevels);


for iRep = 1 : nRep
    %% Compute error 
    [errorTmp] = spaceData.ComputeErrorFunc(bankOfEstimations{iRep}, fucnTest{iRep},...
                                            gridLimitsTest, nDim, false,[], distTest, true);
    errorMean(iRep,:) = mean(abs(errorTmp),1);
    errorMedian(iRep,:) = median(abs(errorTmp),1);
    errorMax(iRep,:) = max(abs(errorTmp), [], 1);

end

end