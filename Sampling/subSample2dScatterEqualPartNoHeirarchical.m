function [samplingEachStage] = subSample2dScatterEqualPartNoHeirarchical(samplingPlaces, ratesToSubSample, nDim)
disp('Notice that your are using a non-heirrchical method');
%% init
nLevels = length(ratesToSubSample);
n = size(samplingPlaces(:,1),1); % number of data points

samplingEachStage = false(n, length(ratesToSubSample));

basicVector = zeros(n,1);

%% 
lastInd = 1;
for index = 1 : nLevels
    basicVector(:) = 0;
    if index == nLevels
        newLastInt = n;
    else
        newLastInt = lastInd + ceil((n / nLevels));
    end
    indces = lastInd : newLastInt; % Notice: I made change here!!
    lastInd = newLastInt;
    
    basicVector(indces) = true;
    samplingEachStage(:, index) = basicVector;   
end
end