function [funcNoised] = AddOutliersAnat(func, pOut, fillDistance, scenario)

[funcNoised] = @(X) SunBfunc(X, func, pOut, fillDistance, scenario);
end

function [values, ps] = SunBfunc(x, func, pOut, h, scenario)
% error('Need to check');
if scenario == 1
    min1 = h^2/2;
    max1 = 3*h^2/2;
elseif scenario == 2
    min1 = h/2;
    max1 = 3*h/2;
elseif scenario == 3
    min1 = 1;
    max1 = 1.5;
end

loc = rand(size(x,1), 1);
ps = loc < pOut; % where outlier should be

neg = rand(size(x,1), 1);
signs = ones(size(x,1), 1);
signs(neg<0.5) = -1;

values =  func(x) + ...
            signs .* ((max1 - min1).* rand(size(x,1), 1) + min1) .* ps;

end