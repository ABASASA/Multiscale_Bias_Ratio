function [errorSource, pointsTest] = ComputeError(fEstimate, func, gridLimits,...
        nDim, flag, outliersPos, testDist, flagAllLevels, flagRelError)

if ~exist('flag', 'var')
    flag = false;
end

if ~exist('outliersPos','var')
    outliersPos = [];
end
  
if ~exist('testDist','var')
    testDist = 0.01;
end

if ~exist('flagAllLevels','var')
    flagAllLevels = true;
end

if ~exist('flagRelError','var')
    flagRelError = 0;
end

if flagRelError == 0 % use relative error
    errorFunc = @(x, y) abs(x - y) ./ x;
elseif flagRelError == 1% use aboulte difference
    errorFunc = @(x, y) abs(x - y);
elseif flagRelError == 2
    errorFunc = @(x, y) x - y;
else
    error('None defined error function');
end

if flag
    figure;
end

if nDim == 2
    [errorSource, pointsTest] = ComputeError2Dim(fEstimate, func, gridLimits, flag,...
    outliersPos, testDist, flagAllLevels, errorFunc, nDim);
elseif nDim == 1
    [errorSource, pointsTest] = ComputeError1Dim(fEstimate, func, gridLimits, flag,...
        outliersPos, testDist, flagAllLevels, errorFunc, nDim);
else
    error('TBD error function');
end
 
end

function [errorSource, pointsTest] = ComputeError2Dim(fEstimate, func, gridLimits, flag,...
    outliersPos, testDist, flagAllLevels, errorFunc, nDim)

% [pointsTest, ~] = Sample2DfunctionEqui(func, gridLimits, testDist);
% [pointsTest, ~] = SampleNDfunctionScatterUniform(func, gridLimits, testDist, nDim);
[pointsTest, ~,~] = SampleRnfunctionEqui(func, gridLimits, testDist, nDim);


% errorSource = zeros(size(pointsTest,1), length(fEstimate)); 
errorScales = zeros(size(pointsTest,1), length(fEstimate));

tmp = zeros(size(pointsTest,1), 1);
sizeSide = sqrt(size(pointsTest,1));
counter = 1;
tmpOld = [];

if flagAllLevels
    startLevel = 1;
    errorSource = zeros(size(pointsTest,1), length(fEstimate)); 
else
    startLevel = length(fEstimate);
    errorSource = zeros(size(pointsTest,1), 1); 
end

for iLevel = startLevel : length(fEstimate)
    tmpFucntion = fEstimate{iLevel};
    
    parfor iPoint = 1 : size(pointsTest,1)
        tmp(iPoint) = tmpFucntion(pointsTest(iPoint,1:2));
    end
    if flagAllLevels

        errorSource(:, iLevel) = errorFunc(pointsTest(:,3), tmp(:));
    else
        errorSource = errorFunc(pointsTest(:,3), tmp(:));
    end
    % Ploting
    if flag
        subplot(length(fEstimate), 3, counter);
        imagesc(pointsTest(:,1), pointsTest(:,2), reshape(tmp, [sizeSide,sizeSide]));
        title(['Level ' num2str(iLevel) ' - result'])
        hold on;
        scatter(outliersPos(:,1), outliersPos(:,2),30,'r');
        colorbar;

        subplot(length(fEstimate), 3, counter + 1);
        if flagAllLevels
            imagesc(pointsTest(:,1), pointsTest(:,2),reshape(errorSource(:, iLevel), [sizeSide,sizeSide]));
        else
            imagesc(pointsTest(:,1), pointsTest(:,2),reshape(errorSource(:), [sizeSide,sizeSide]));
        end
        hold on;
        scatter(outliersPos(:,1), outliersPos(:,2),30,'r');
%         plot3(pointsTest(:,1), pointsTest(:,2), errorSource(:, iLevel))

        title(['Level ' num2str(iLevel) ' - Source. Rell err_{mean} = ' ...
                        num2str(mean((errorSource(:, iLevel))))])
        colorbar;

    end
    
    if (iLevel ~= 1 && flagAllLevels)
        errorScales(:, iLevel - 1) = errorFunc(oldTmp(:), tmp(:));
        % Ploting
        if flag
             subplot(length(fEstimate), 3, counter + 2);
            imagesc(pointsTest(:,1), pointsTest(:,2),...
                reshape(errorScales(:, iLevel - 1), [sizeSide,sizeSide]));
            hold on;
            scatter(outliersPos(:,1), outliersPos(:,2),30,'r');
            title(['Level ' num2str(iLevel) ' - Scales . Rel err_{mean} = ' ...
                            num2str(mean((errorScales(:, iLevel-1))))])
            colorbar;
            
            if iLevel == length(fEstimate)
                figure;
                imagesc(pointsTest(:,1), pointsTest(:,2),...
                        reshape(errorScales(:, iLevel - 1), [sizeSide,sizeSide]));
                hold on;
                scatter(outliersPos(:,1), outliersPos(:,2),30,'r');
                title(['Level ' num2str(iLevel) ' - Scales . Rel err_{mean} = ' ...
                                num2str(mean((errorScales(:, iLevel-1))))])
                colorbar;
                figure;
                imagesc(pointsTest(:,1), pointsTest(:,2), reshape(tmp, [sizeSide,sizeSide]));
                title(['result'])
                hold on;
                scatter(outliersPos(:,1), outliersPos(:,2),30,'r');
                colorbar;
            end
            
        end
        
        
    end
    oldTmp = tmp;
    
    counter = counter + 3;

end

end


function [errorSource, pointsTest] = ComputeError1Dim(fEstimate, func, gridLimits, flag,...
    outliersPos, testDist, flagAllLevels, errorFunc, nDim)

[pointsTest, ~] = Sample1DfunctionEqui(func, gridLimits, testDist);
% errorSource = zeros(size(pointsTest,1), length(fEstimate)); 
errorScales = zeros(size(pointsTest,1), length(fEstimate));

tmp = zeros(size(pointsTest,1), 1);
sizeSide = sqrt(size(pointsTest,1));
counter = 1;
tmpOld = [];

if flagAllLevels
    startLevel = 1;
    errorSource = zeros(size(pointsTest,1), length(fEstimate)); 
else
    startLevel = length(fEstimate);
    errorSource = zeros(size(pointsTest,1), 1); 
end

for iLevel = startLevel : length(fEstimate)
    tmpFucntion = fEstimate{iLevel};
    
    parfor iPoint = 1 : size(pointsTest,1)
        tmp(iPoint) = tmpFucntion(pointsTest(iPoint,1:nDim));
    end
    if flagAllLevels

        errorSource(:, iLevel) = errorFunc(pointsTest(:,nDim+1), tmp(:));
    else
        errorSource = errorFunc(pointsTest(:,nDim+1), tmp(:));
    end
    % Ploting
    if flag
        subplot(length(fEstimate), 3, counter);
        plot(pointsTest(:,1), tmp, '*-');
        title(['Level ' num2str(iLevel) ' - result'])
        hold on;
        scatter(outliersPos(:,1), outliersPos(:,2),30,'r');

        subplot(length(fEstimate), 3, counter + 1);
        if flagAllLevels
            plot(pointsTest(:,1), errorSource(:, iLevel), '*-');
        else
            plot(pointsTest(:,1), errorSource(:), '*-');
        end
        hold on;
        scatter(outliersPos(:,1), outliersPos(:,2),30,'r');
%         plot3(pointsTest(:,1), pointsTest(:,2), errorSource(:, iLevel))

        title(['Level ' num2str(iLevel) ' - Source. Rell err_{mean} = ' ...
                        num2str(mean((errorSource(:, iLevel))))])

    end
    
    if (iLevel ~= 1 && flagAllLevels)
        errorScales(:, iLevel - 1) = errorFunc(oldTmp(:), tmp(:));
        % Ploting
        if flag
             subplot(length(fEstimate), 3, counter + 2);
              plot(pointsTest(:,1), errorScales(:, iLevel - 1), '*-');
            hold on;
            scatter(outliersPos(:,1), outliersPos(:,2),30,'r');
            title(['Level ' num2str(iLevel) ' - Scales . Rel err_{mean} = ' ...
                            num2str(mean((errorScales(:, iLevel-1))))])
            
            if iLevel == length(fEstimate)
                figure;
                plot(pointsTest(:,1), errorScales(:, iLevel - 1), '*-');
                hold on;
                scatter(outliersPos(:,1), outliersPos(:,2),30,'r');
                title(['Level ' num2str(iLevel) ' - Scales . Rel err_{mean} = ' ...
                                num2str(mean((errorScales(:, iLevel-1))))])
                figure;
                plot(pointsTest(:,1), tmp, '*-');
                title(['result'])
                hold on;
                scatter(outliersPos(:,1), outliersPos(:,2),30,'r');
            end
            
        end
        
        
    end
    oldTmp = tmp;
    
    counter = counter + 3;

end

end

