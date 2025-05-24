function [errorMat, PointTest] =...
        ComputeErrorAccording3MeasuremntsByPoint(bankOfEstimations, fucnTest,...
                                          gridLimitsTest, distTest, nRep, nDim)


nLevels = length(bankOfEstimations{1});


for iLevel = 1 : nLevels
    funcBank = cell(nRep,1);
    for iRep = 1 : nRep
        tmpp = bankOfEstimations{iRep};
        funcBank{iRep} = tmpp{iLevel};
    end
    %% Compute error 
    [errorTmp, PointTest] = ComputeError(funcBank, fucnTest,...
                gridLimitsTest, nDim, false,[], distTest, true);
%     if nDim == 1
        if iLevel == 1
%             errorMean = zeros(size(PointTest,1), nLevels);
%             errorMedian = zeros(size(PointTest,1), nLevels);
%             errorMax = zeros(size(PointTest,1), nLevels);
            errorMat = zeros(size(PointTest,1), nLevels, 3);
        end
%     else
%         error('1D only');
%     end

    errorMat(:,iLevel, 1) = mean(abs(errorTmp),2);
    errorMat(:,iLevel, 2) = median(abs(errorTmp),2);
    errorMat(:,iLevel, 3) = max(abs(errorTmp), [], 2);

end
PointTest = PointTest(:,1);

end