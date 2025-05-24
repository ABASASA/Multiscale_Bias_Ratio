function [func] = CombineNonSmoothRPlus(g,r)
rr = @(x,y) RPlus(r, x, y);
func = @(x,y) g(x,y) + rr(x,y);


end

function [rp] = RPlus(r, x, y)
rp = 0 * x;

v = r(x,y);

rp(v >=0) = v(v>=0);


end