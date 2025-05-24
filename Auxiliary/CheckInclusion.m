function [] = CheckInclusion(samplingEachStage)

nLevels = size(samplingEachStage, 2);
flags = false(nLevels, 1);
for iLevel = 1 : nLevels - 1
    
    A = samplingEachStage(:, iLevel+ 1) - samplingEachStage(:, iLevel);
    if sum(A < 0)
        flags(iLevel) = true;
    end
    
end

inds = 1 : nLevels;
inds = inds(flags);

disp(['Non-inclusion in levels: ' num2str(inds)]);


end