function [samplingEachStage] = subSample2dScatterWithInclusion(samplingPlaces, ratesToSubSample)

% init
n = length(samplingPlaces);
samplingEachStage = false(n, length(ratesToSubSample));

basicVector = samplingPlaces(:,1);
basicVector(:) = 0;

for index = 1 : length(ratesToSubSample)
    if index == 1
        basicVector(:) = 0;
    else
        basicVector = samplingEachStage(:, index - 1);
    end
    
    % The numer of new inces to sample
    nIndToSample = round(sqrt(n) / ratesToSubSample(index)).^2 - sum(basicVector);
    
    
    indces = randperm(n - sum(basicVector), nIndToSample);
    subBasicVector = zeros(n - sum(basicVector), 1);
    subBasicVector(indces) = true;
    basicVector(~basicVector) = subBasicVector;
    
    
%     if sum(indces == size(basicVector,2)) == 0
%         indces(end + 1) = length(basicVector);
%     end
%     basicVector(indces) = true;
    samplingEachStage(:, index) = basicVector;
    
end
end