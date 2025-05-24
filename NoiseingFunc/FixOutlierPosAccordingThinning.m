function [outliersPosQueued] = FixOutlierPosAccordingThinning(outliersPos, queduedIndex)
N = length(outliersPos);
indces = 1 : N;
outIndeces = indces(outliersPos);

outliersPosQueued = false(N,1); % init
for i = 1 : length(outIndeces)
    outliersPosQueued(queduedIndex == outIndeces(i)) = true;
end

end