% datPoint is the total set of Points currently assum it is in [-10,10]^2:
% (x, y, value).
function [f,e, additionalInput] = MultiScale(dataPoints,samplingEachStage,...
                                            kernel, deltas, estFunc, nDim, spaceData)
                                        
[f,e, additionalInput] = MultiScaleV1(dataPoints,samplingEachStage,...
                                            kernel, deltas, estFunc, nDim, spaceData);
end

function [f,e, additionalInput] = MultiScaleV1(dataPoints,samplingEachStage,...
                                            kernel, deltas, estFunc, nDim, spaceData)

minusFunc = spaceData.minusFunc;
plusFunc = spaceData.plusFunc;
baseElement = spaceData.baseElement;

N = length(deltas);

e = cell(N,1);
f = cell(N,1);
s = cell(N,1);
additionalInput = cell(N,1);

for iLevel = 1 : N
    %% Current set of points
    delta = deltas(iLevel);
    Xj = dataPoints(samplingEachStage(:,iLevel),:);
    
    %% Compute projection error
    if iLevel == 1
        e{iLevel} = minusFunc(Xj(:, nDim + 1 : end), baseElement);
%         e{iLevel} = Xj(:, nDim + 1 : end);

    else
%         reminder = Xj(:, nDim + 1 : end);
%         for k = 1 : iLevel-1
%             reminder = reminder - LoopOverSFunction(s{k},Xj(:,1:nDim), size(Xj,2) - nDim);
%         end  
        e{iLevel} = minusFunc(Xj(:, nDim + 1 : end), LoopOverSFunction(f{iLevel - 1}, ...
                                    Xj(:,1:nDim), size(Xj,2) - nDim));
    end

    XjMod = Xj;
    XjMod(:, nDim + 1 : end) = plusFunc(e{iLevel}, baseElement) ; % Update the values to reminder
    e{iLevel};
    %% estimation
    [s{iLevel}, additionalInput{iLevel}] = estFunc(XjMod, delta, kernel,...
                                            additionalInput, nDim);   
    %% Compute f
    if (iLevel == 1)
        f{iLevel} = @(X) plusFunc(minusFunc(s{iLevel}(X), baseElement), baseElement);
    else
        f{iLevel} = @(X) plusFunc(minusFunc(s{iLevel}(X), baseElement), f{iLevel-1}(X));
    end

    %% disp
    disp(['Finish Level: ' num2str(iLevel)]);
end

end






