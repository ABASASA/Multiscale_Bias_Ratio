function [spaceData] = SpaceDataSPD3()
spaceData = struct;

spaceData.minusFunc = @(Xs, Y) LogSPD(Xs, Y);
spaceData.plusFunc = @(Xs, Y) ExpSPD(Xs, Y);
spaceData.baseElement = eye(3,3);
spaceData.ComputeErrorFunc = @ComputeErrorSPD;
spaceData.geoDist = @GeodesicDistanceSPD;
end
