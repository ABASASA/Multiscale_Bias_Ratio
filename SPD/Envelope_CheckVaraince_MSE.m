function [spaceData] = SpaceDataSPD3()
    spaceData = struct;
    spaceData.minusFunc = @(Xs, Y) LogSPD(Xs, Y);
    spaceData.plusFunc = @(Xs, Y) EXPSPD(Xs, Y);
    spaceData.baseElement = eye(3,3);
end
                                     