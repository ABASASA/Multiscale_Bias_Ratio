function [samplingEachStage] = subSample2dRandom(samplingPlaces, ratesToSubSample)

% init
N = size(samplingPlaces,1);

samplingEachStage = false(N, length(ratesToSubSample));

basicVector = false(N,1);

randVec = rand(N, 1);

%% Main code.
% The idea is to use the same datapoint in previus levels with an additonal
% points. We intialize a uniform vector betweem [0,1] in length N. The
% 1/rate is the threshold on values which lower than it we will use and
% abouve do not. E.g, for 1/1 we use the whole data. The data size is
% increases as the rate increase.

for index = 1 : length(ratesToSubSample)
    
    percentOfData = (1 / ratesToSubSample(index));
    basicVector(randVec <= percentOfData) = true;
    samplingEachStage(:, index) =  basicVector;
    
end
end