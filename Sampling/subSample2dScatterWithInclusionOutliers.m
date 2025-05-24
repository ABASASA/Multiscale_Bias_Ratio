function [samplingEachStage] = subSample2dScatterWithInclusionOutliers(...
    samplingPlaces, ratesToSubSample, outliersPos)

% init
n = length(samplingPlaces);
nOutliers = sum(outliersPos);

samplingEachStage = false(n, length(ratesToSubSample));

basicVector = samplingPlaces(~outliersPos,1);
basicVector(:) = 0;

for index = 1 : length(ratesToSubSample)-1
    if index == 1
        basicVector(:) = 0;
    else
        basicVector = samplingEachStage(~outliersPos, index - 1);
    end
    
    % The numer of new inces to sample
    nIndToSample = round(sqrt(n - nOutliers) / ratesToSubSample(index)).^2 - sum(basicVector);
    
    length(basicVector) - sum(basicVector) 
    nIndToSample
    indces = randperm(length(basicVector) - sum(basicVector), nIndToSample);
    subBasicVector = zeros(length(basicVector) - sum(basicVector), 1);
    subBasicVector(indces) = true;
    basicVector(~basicVector) = subBasicVector;
    
    
%     if sum(indces == size(basicVector,2)) == 0
%         indces(end + 1) = length(basicVector);
%     end
%     basicVector(indces) = true;
    samplingEachStage(~outliersPos, index) = basicVector;
    
end
samplingEachStage(:, end) = true;
end