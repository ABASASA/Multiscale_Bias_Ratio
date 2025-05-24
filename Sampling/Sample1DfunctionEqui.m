function [pointsSet, samplingPlaces, outliersPos] = Sample1DfunctionEqui(...
                func, sizes, distAxis)

samplingPlaces = sizes(1) :  distAxis : sizes(2);
samplingPlaces = samplingPlaces';

[value, outliersPos] = func(samplingPlaces);
% pointsSet = [samplingPlaces, value];

pointsSet = [samplingPlaces, value];

end