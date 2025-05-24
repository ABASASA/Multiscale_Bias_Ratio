function [pointsSet, samplingPlaces, outliersPos] = SampleRnfunctionEqui(...
                func, sizes, distAxis, nDim)
if nDim == 2
    [pointsSet, samplingPlaces, outliersPos]  = SampleR2functionEqui(func,...
    sizes, distAxis);
elseif nDim == 1
    [pointsSet, samplingPlaces, outliersPos]  = SampleR1functionEqui(func,...
    sizes, distAxis);
else
    error('Only for 1 and 2 dim');
end
end

function [pointsSet, samplingPlaces, outliersPos]  = SampleR1functionEqui(func,...
    sizes, distAxis)
samplingPlaces = sizes(1) :  distAxis : sizes(2);
samplingPlaces = transpose(samplingPlaces);

[value, outliersPos] = func(samplingPlaces);
% pointsSet = [samplingPlaces, value];

pointsSet = [samplingPlaces, value];
end

function [pointsSet, samplingPlaces, outliersPos]  = SampleR2functionEqui(func,...
    sizes, distAxis)
samplingPlaces = sizes(1) :  distAxis : sizes(2);
[X,Y] = meshgrid(samplingPlaces, samplingPlaces);
X = X(:);
Y = Y(:);

[value, outliersPos] = func([X, Y]);
% pointsSet = [samplingPlaces, value];

pointsSet = [X, Y, value];
samplingPlaces = [X, Y];
end