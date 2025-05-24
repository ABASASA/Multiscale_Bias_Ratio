function [error] = ComputeErrorBeweenScales(fEstimate, func, gridLimits, flag)

if ~exist('flag', 'var')
    flag = false;
end
nScales = length(fEstimate);
if nScales <= 1
    disp('Not enough scales')
    return;
end

[pointsTest, samplingPlaces] = Sample2DfunctionEqui(func, gridLimits, 0.005);
error = zeros(size(pointsTest,1), length(fEstimate) - 1);
tmp = zeros(size(pointsTest,1), 1);
sizeSide = sqrt(size(pointsTest,1));
counter = 1;
oldTmp = [];


figure;

for iLevel = 1 : length(fEstimate)
    tmpFucntion = fEstimate{iLevel};
    
    parfor iPoint = 1 : size(pointsTest,1)
        tmp(iPoint) = tmpFucntion(pointsTest(iPoint,1),pointsTest(iPoint,2));
    end
    
    if flag
        subplot(length(fEstimate), 2, counter);
        imagesc(reshape(tmp, [sizeSide,sizeSide]));
        title(['Level ' num2str(iLevel) ' - result'])
        colorbar;

    end
    

    if (iLevel ~= 1)
        error(:, iLevel - 1) = (oldTmp(:) - tmp(:));

        % Ploting
        if flag
             subplot(length(fEstimate), 2, counter + 1);
            imagesc(reshape(error(:, iLevel - 1), [sizeSide,sizeSide]));
            title(['Level ' num2str(iLevel) ' - diffrence. err_{mean} = ' ...
                            num2str(mean(abs(error(:, iLevel-1))))])
            colorbar;

        end
    end
    oldTmp = tmp;
    counter = counter + 2;

        
end
end