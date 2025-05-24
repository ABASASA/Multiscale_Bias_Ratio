function [funcNoised] = AddGuassianNoiseAndOutliers(func, sigma, pOut, sigmaOut)

[funcNoised] = @(x) SunBfunc(x, func, sigma, pOut, sigmaOut);
end

function [values, ps] = SunBfunc(x, func, sigma, pOut, sigmaOut)
% The noise is proportional to the function value
ps = rand(size(x,1),1) < pOut;
fx = func(x);
values =  fx + (randn(size(x,1), 1) .* sigma .* fx)  + ...
                (fx .* sigmaOut .* randn(size(x,1), 1) .*ps);

end

function [values, ps] = SunBfuncAbsoulte(x, func, sigma, pOut, sigmaOut)
% The noise is absoulte
ps = rand(size(x,1),1) < pOut;
values =  func(x) + (randn(size(x,1), 1) .* sigma)  + ...
                (func(x) .* sigmaOut .* randn(size(x,1), 1) .*ps);

end