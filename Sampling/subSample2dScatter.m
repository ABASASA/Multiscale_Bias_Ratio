function [samplingEachStage] = subSample2dScatter(samplingPlaces, ratesToSubSample)

% init
n = length(samplingPlaces);
samplingEachStage = false(n, length(ratesToSubSample));

basicVector = samplingPlaces(:,1);
basicVector(:) = 0;

for index = 1 : length(ratesToSubSample)
    basicVector(:) = 0;

    indces = randperm(n, round(sqrt(n) / ratesToSubSample(index)).^2);
    
    
%     if sum(indces == size(basicVector,2)) == 0
%         indces(end + 1) = length(basicVector);
%     end
    basicVector(indces) = true;
    samplingEachStage(:, index) = basicVector;
    
end
end