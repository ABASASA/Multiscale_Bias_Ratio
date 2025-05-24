function [pointsSet, samplingPlaces, outliersPos] = SampleNDfunctionScatterUniform(func,...
                                                        sizes, distAxis, nDim)

nSamples = length(sizes(1) :  distAxis : sizes(2)).^nDim; % equal number to equidistance

samplingPlaces = sizes(1) + (sizes(2) - sizes(1)) * rand(nSamples, nDim);

% [value, outliersPos] = func(samplingPlaces(:, 1), samplingPlaces(:, 2));
[value, outliersPos] = func(samplingPlaces);

pointsSet = [samplingPlaces, value];

end