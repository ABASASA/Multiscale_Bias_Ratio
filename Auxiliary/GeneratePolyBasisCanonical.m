function [poly] = GeneratePolyBasisCanonical(m)

ind = 0 : m;
n = length(ind);
XSpace = repmat(ind, 1, n);
XSpace = XSpace(:);
YSpace = repmat(ind, n, 1);
YSpace = YSpace(:);

poly = @(x,y) transpose((x.^XSpace) .* (y.^YSpace));

end

