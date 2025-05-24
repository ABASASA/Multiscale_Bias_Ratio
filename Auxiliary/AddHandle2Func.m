function [funcHand] = AddHandle2Func(func)
% This function created to replace 'AddGuassianNoiseAndOutliers' in the places where
% noise is not needed.
    funcHand = @(X) TmpFunc(func, X);
end

function [values, outliersPos] = TmpFunc(func, X)
    values = func(X);
    outliersPos = false(size(X,1),1);
end