function [queduedPointSet, samplingPlacesQueued, outliersPosQueued] = Sample2DfunctionScatterUniformWithThinningGreedy(func,...
                                                        sizes, distAxis)

nSamples = length(sizes(1) :  distAxis : sizes(2)).^2; % equal number to equidistance

samplingPlaces = sizes(1) + (sizes(2) - sizes(1)) * rand(nSamples,2);

[value, outliersPos] = func(samplingPlaces(:, 1), samplingPlaces(:, 2));

pointsSet = [samplingPlaces, value];

[queduedPointSet, queduedIndex, ~] = GreedyThiningBroneWende(pointsSet);
[outliersPosQueued] = FixOutlierPosAccordingThinning(outliersPos, queduedIndex);
samplingPlacesQueued = queduedPointSet(:,1:2);
end
