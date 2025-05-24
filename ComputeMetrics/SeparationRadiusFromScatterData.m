function [distances] = SeparationRadiusFromScatterData(pointsSet, samplingEachStage)

nLevels = size(samplingEachStage, 2);
distances = zeros(nLevels, 1);

for  i = 1 : nLevels
   Xj = pointsSet(samplingEachStage(:,i), :);
   [distances(i), ~] = ComputeSeparationRadiusRn(Xj);

end

end