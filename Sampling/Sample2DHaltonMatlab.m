function [pointsSet, samplingPlaces, outliersPos] = Sample2DHaltonMatlab(func,...
                                                        gridLimits, distAxis)

p = haltonset(2,'leap',1000,'skip',100);                                                    
% P is from [0,1] and we only allow between [0, 1 - dist]

nSamples = length(gridLimits(1) :  distAxis : gridLimits(2)).^2; % equal number to equidistance
% samplesInd = randperm(size(p,1), nSamples);
samplesInd = 1 : nSamples;
samplingPlaces = p(samplesInd, :);
samplingPlaces = (gridLimits(2) - gridLimits(1)) .* samplingPlaces + gridLimits(1);

[value, outliersPos] = func(samplingPlaces(:, 1), samplingPlaces(:, 2));
pointsSet = [samplingPlaces, value];

end