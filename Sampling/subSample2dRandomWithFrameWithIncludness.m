function [samplingEachStage] = subSample2dRandomWithFrameWithIncludness(samplingPlaces, ratesToSubSample)

% init
samplingEachStage = false(length(samplingPlaces).^2, length(ratesToSubSample));

basicVectorEquidistnace = samplingPlaces;


for index = length(ratesToSubSample) : -1 : 1
    basicVectorEquidistnace(:) = 0;

    indces = zeros(length(1 : ratesToSubSample(index) : length(basicVectorEquidistnace)) +1 , 1);
    indces(1) = 1;
    indces(end) = length(basicVectorEquidistnace);
    
    if index == length(ratesToSubSample)
        indces(2:end-1) = randi([2, length(basicVectorEquidistnace)-1],...
                                                [length(indces) - 2,1] );
        maxSet = indces;
    else
        tmp = randi([2, length(maxSet)-1], [length(indces) - 2,1] );
        indces(2:end-1) = maxSet(tmp);
        maxSet = indces;
    end
    basicVectorEquidistnace(indces) = true;

    [X,Y] = meshgrid(basicVectorEquidistnace, basicVectorEquidistnace);
    X = X .* Y;
    samplingEachStage(:, index) =  X(:);
    
end
end