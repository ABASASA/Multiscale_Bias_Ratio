function [errorMean, errorMedian, errorMax ] =...
        ComputeErrorAccording3Measuremnts(bankOfEstimations, fucnTest,...
                                          gridLimitsTest, distTest, nRep, nDim)

errorMean = zeros(nRep, 1);
errorMedian = zeros(nRep, 1);
errorMax = zeros(nRep, 1);

for iRep = 1 : nRep
    %% Compute error 
    [errorTmp] = ComputeError(bankOfEstimations{iRep}, fucnTest{iRep},...
                                            gridLimitsTest, nDim, false,[], distTest, false);
    errorMean(iRep) = mean(abs(errorTmp));
    errorMedian(iRep) = median(abs(errorTmp));
    errorMax(iRep) = max(abs(errorTmp));
end

end