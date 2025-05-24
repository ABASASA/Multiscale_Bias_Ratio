function [samplingEachStage] = subSample2dRandomWithFrame(samplingPlaces, ratesToSubSample)

% init
samplingEachStage = false(length(samplingPlaces).^2, length(ratesToSubSample));

basicVectorEquidistnace = samplingPlaces;
basicVectorEquidistnace(:) = 0;


for index = 1 : length(ratesToSubSample)
    basicVectorEquidistnace(:) = 0;
    indces = zeros(length(1 : ratesToSubSample(index) : length(basicVectorEquidistnace)) +1 , 1);
    indces(1) = 1;
    indces(end) = length(basicVectorEquidistnace);
    
    
    indces(2:end-1) = randi([2, length(basicVectorEquidistnace)-1], [length(indces) - 2,1] );
    basicVectorEquidistnace(indces) = true;

    [X,Y] = meshgrid(basicVectorEquidistnace, basicVectorEquidistnace);
    X = X .* Y;
    samplingEachStage(:, index) =  X(:);
    
end
end