function [norms] = MeshNormFromScatterData(pointsSet, samplingEachStage, gridLimits)

nLevels = size(samplingEachStage, 2);
norms = zeros(nLevels, 1);

for  i = 1 : nLevels
   Xj = pointsSet(samplingEachStage(:,i), :);
   norms(i) = ComputeMeshNorm(Xj(:,1), Xj(:,2), gridLimits);
end

end