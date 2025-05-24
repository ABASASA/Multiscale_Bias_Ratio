function [spaceData] = SpaceDataEuclidian()

spaceData = struct;

spaceData.minusFunc = @(Xs, Y) Xs - Y;
spaceData.plusFunc = @(Xs, Y) Xs + Y;
spaceData.baseElement = 0;
spaceData.ComputeErrorFunc = @ComputeError;

end