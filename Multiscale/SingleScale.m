function [f,e, additionalInput] = SingleScale(dataPoints, ...
                        samplingEachStage, kernel, deltas, estFunc, nDim, spaceData)

N = length(deltas);

e = cell(N,1);
f = cell(N,1);

additionalInput = cell(N,1);
%% In each number of data chooce the relevent one.
for iLevel = 1 : N
    delta = deltas(iLevel);
    samplingEachStage1 = samplingEachStage(:, iLevel);

    [f{iLevel},e{iLevel}, additionalInput{iLevel}] = MultiScale(dataPoints,...
                        samplingEachStage1, kernel, delta, estFunc, nDim, spaceData);
    f{iLevel} = f{iLevel}{1};
    e{iLevel} = e{iLevel}{1};
    additionalInput{iLevel} = additionalInput{iLevel}{1};
    
end