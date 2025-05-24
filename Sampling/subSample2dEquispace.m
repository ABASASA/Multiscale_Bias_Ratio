function [samplingEachStage] = subSample2dEquispace(samplingPlaces, ratesToSubSample)

% init
samplingEachStage = false(length(samplingPlaces), length(ratesToSubSample));

basicVector = samplingPlaces;
basicVector(:) = 0;


for index = 1 : length(ratesToSubSample)
    indces = 1 : ratesToSubSample(index) : length(basicVector);
%     if sum(indces == size(basicVector,2)) == 0
%         indces(end + 1) = length(basicVector);
%     end
    basicVector(indces) = true;
    [X,Y] = meshgrid(basicVector, basicVector);
    X = X .* Y;
    samplingEachStage(:, index) =  X(:);
    
end
end