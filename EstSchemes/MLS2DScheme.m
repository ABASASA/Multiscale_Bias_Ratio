function [currentEstimation, additionalInput] = MLS2DScheme(XjMod, delta,...
                    kernel, additionalInput, nDim, m) 
if nDim == 2
    currentEstimation = @(x) MLS2D(XjMod, m, kernel, delta, x);
else
    erorr('I did not implement this dimension');
end

end