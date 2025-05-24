function [measurementPerPoint, measurementIntegrated] = ComputeMSEMeusurmentesSPD(...
            bankOfEstimations, funcTestforBias, testPointss,...
            nRep, nDimOutPut, spaceData, intergrationFactor)

    nLevels = length(bankOfEstimations{1}); % number of levels
    nTestPointss = size(testPointss, 1); % number of points in the test grid
    measurementPerPoint = zeros(nLevels, nTestPointss, 3);
    measurementIntegrated = zeros(nLevels, 3);
    
    %% Compute each compnenet in each level
    for iLevel = 1 : nLevels
        
        %% First - We build a list of all relvent functions
        listFunc = cell(nRep,1);
        fMean = @(x) 0;
        for jFunc = 1 : nRep % compute the error for each function then mean
            func = GetFunc(bankOfEstimations, jFunc, iLevel);
            listFunc{jFunc} = func;
            if iscell(func)
                error('Problem Cell');
            end
        end
        fMean = @(point) ComputeMean(point, listFunc, nDimOutPut, spaceData);

        %% Compute the MSE and its properties in each level 
        [measurementPerPoint(iLevel, :, 1), measurementPerPoint(iLevel, :, 2), ...
            measurementPerPoint(iLevel, :, 3), measurementIntegrated(iLevel, 1), ...
            measurementIntegrated(iLevel, 2), measurementIntegrated(iLevel, 3)]  =...
                    ComputeErrorMSE(listFunc, fMean, funcTestforBias,...
                                        testPointss, intergrationFactor, spaceData);
         
    end
        
end

function [MISE_I, Bias_I, Variance_I, IMSE, IBiasSquared,...
        IVar] = ComputeErrorMSE(funcList, fMean, funcTestforBias, ...
                            testPointss, intergrationFactor, spaceData)

    nTestPoints = size(testPointss, 1);
    nRep = length(funcList);
    MISE_I = zeros(nTestPoints, 1);
    Bias_I = zeros(nTestPoints, 1);
    Variance_I = zeros(nTestPoints, 1);
    midCompIVar = zeros(nTestPoints, nRep);
    midCompIMSE = zeros(nTestPoints, nRep);
    % For each point in the test set it compute the variance, MSE and bias
    parfor i = 1 : nTestPoints
        
        point = testPointss(i,:); % Current Point
       
        fMean_I = fMean(point); % FMean value
        base = funcTestforBias(point); % ground truth
        
        tmpMSE = zeros(nRep, 1);
        tmpVariance = zeros(nRep, 1);

        for j = 1 : nRep % Go over every rep in the monte carelo
            func = funcList{j};
            funcValue = func(point);
            tmpMSE(j) = spaceData.geoDist(funcValue, base) .^ 2;
            tmpVariance(j) = spaceData.geoDist(funcValue, fMean_I) .^ 2;
        end
        MISE_I(i, :) = mean(tmpMSE);
        
        Bias_I(i) = spaceData.geoDist(fMean_I, funcTestforBias(testPointss(i,:))) .^ 2;
        
        Variance_I(i) = mean(tmpVariance);      
        midCompIVar(i, :) = tmpVariance;
        midCompIMSE(i, :) = tmpMSE;
    end
    %% compute the Integreted measuerment
    IBiasSquared = sum(Bias_I) * intergrationFactor;
    IVar = mean(intergrationFactor * sum(midCompIVar, 1));
    IMSE = mean(intergrationFactor * sum(midCompIMSE, 1));
end


function [means] = ComputeMean(point, listFunc, nDimOutPut, spaceData)

%% Hard coded input
nIter = 15;

%% initizalize
nRep = length(listFunc);
values = zeros(nRep, nDimOutPut);

% Precomputing the valuse.
for ii = 1 : nRep
    funcTmp = listFunc{ii};
    values(ii, :) = funcTmp(point);
end
%% Compute Mean

currentGuess = values(1, :);

logFunc = spaceData.minusFunc;
expFunc = spaceData.plusFunc;

for ii = 2 : nIter + 1
   logs = sum(logFunc(values, currentGuess), 1) ./ nRep;
   currentGuess = expFunc(logs, currentGuess);
end

means = currentGuess;
end
