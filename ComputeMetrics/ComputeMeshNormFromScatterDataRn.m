function [norms] = ComputeMeshNormFromScatterDataRn(samplingPlaces,...
                                    samplingEachStage, gridLimits, nDim)

nLevels = size(samplingEachStage, 2);
norms = zeros(nLevels, 1);
for  i = 1 : nLevels
   Xj = samplingPlaces(samplingEachStage(:,i), :);
   norms(i) = ComputeMeshNormRn(Xj, gridLimits, nDim);

end

end