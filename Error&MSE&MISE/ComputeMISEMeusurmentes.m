function [MISE, biasMISE, vairnaceMISE, variance_Lambda_MISE] = ComputeMISEMeusurmentes(...
            bankOfEstimations, funcTestforBias, nRep, nDim,...
                    gridSize, lambda, flagRHKS, valueNormk)
    nSamples = 1000;%number of function in Monte carlo integration
    nLevels = length(bankOfEstimations{1}); % number of levels
    
    MISE = zeros(nLevels,1);
    biasMISE = zeros(nLevels,1);
    vairnaceMISE = zeros(nLevels,1);
    variance_Lambda_MISE = nan(nLevels,1);
    %% Compute each compnenet in each level
    for iLevel = 1 : nLevels
        
        %% First - We build a list of all relvent functions
        listFunc = cell(nRep,1);
        fMean = @(x) 0;
        for jFunc = 1 : nRep % compute the error for each function then mean
            func = GetFunc(bankOfEstimations, jFunc, iLevel);
            listFunc{jFunc} = func;
            fMean = @(X) fMean(X) + func(X) ./ nRep;
        end
        
        %% MISE 
        [MISE(iLevel), biasMISE(iLevel), vairnaceMISE(iLevel) , varL2]  =...
                    ComputeErrorMISE(listFunc, fMean, funcTestforBias , nDim,...
                                        gridSize, nSamples);
       
        %% Variance Lambda
        if flagRHKS      
            variance_Lambda_MISE(iLevel) = mean( varL2.' + ...
                    valueNormk(:, iLevel) .* lambda(:, iLevel));
        end
         
          
    end
        
end

function [MISE, Bias, Variance, variance_J] = ComputeErrorMISE(funcList, fMean, funcTestforBias, ...
                            nDim, gridSize, nSamples)
    if ~isvector(gridSize) | length(gridSize) ~=2
        error('User Error: Problem with grid size');
    end
    
    samples = rand(nSamples, nDim) .* (gridSize(2) - gridSize(1)) + gridSize(1);
    nRep = length(funcList);
    MISE_I_J = zeros(nSamples, nRep);
    Bias_I = zeros(nSamples, 1);
    Variance_I_J = zeros(nSamples, nRep);

    parfor i = 1 : nSamples
        
        point = samples(i,:); % Current Point
       
        fMean_I = fMean(point); % FMean value
        base = funcTestforBias(point);
        
        tmpMISE = zeros(nRep, 1);
        tmpVariance = zeros(nRep, 1);

        for j = 1 : nRep
            func = funcList{j};
            funcValue = func(point);
            tmpMISE(j) = (funcValue - base).^2;
            tmpVariance(j) = (funcValue - fMean_I).^2;
        end
        MISE_I_J(i, :) = tmpMISE;
        
        Bias_I(i) = (fMean_I - funcTestforBias(samples(i,:))).^2;
        
        Variance_I_J(i, :) = tmpVariance;      
    end
    
    MISE_J = mean(MISE_I_J,1); % This is an approximation of the integration using Monte CArlo apprxomation for integrals.
    MISE = mean(MISE_J);
    
    Bias = mean(Bias_I);

    variance_J = mean(Variance_I_J,1); % This is an approximation of the integration using Monte CArlo apprxomation for integrals.
    Variance = mean(variance_J);
end
