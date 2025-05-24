function [samplingEachStage] = subSample2dScatterFirst(samplingPlaces, ratesToSubSample, nDim)

%% init
n = size(samplingPlaces(:,1),1); % number of data points

samplingEachStage = false(n, length(ratesToSubSample));

basicVector = zeros(n,1);

%% 
% We assume that the order of the data is random-uniform.
% Hence, in order to sample randomly from it, we need to take the first n
% entries per stage of the multiscale.
for index = 1 : length(ratesToSubSample)
    basicVector(:) = 0;

    indces = 1 : (round(sqrt(n) / ratesToSubSample(index))^nDim); % Notice: I made change here!!
    
    basicVector(indces) = true;
    samplingEachStage(:, index) = basicVector;   
end
end