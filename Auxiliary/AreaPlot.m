function [dataFlip] = AreaPlot(Xdata, curvelower, curveupper)
curve1 = curveupper;
curve2 = curvelower;
X = [Xdata, fliplr(Xdata)];

inBetween = [curve1, fliplr(curve2)];

dataFlip = struct;

dataFlip.X = X;
dataFlip.Y = inBetween;

end