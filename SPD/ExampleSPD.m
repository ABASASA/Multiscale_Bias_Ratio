function [value] = ExampleSPD(X)

A =  @(x,y) [sin(2*pi*y) .* cos(2*pi*x), y.^2, x.*y;  y.^2, 1, 0;  x.*y, 0,  cos(pi * x)];
Nel = size(X,1);
value = zeros(Nel, 9);

for i = 1 : Nel
   tmpp = expm(A(X(i,1),X(i,2)));
   value(i,:) = tmpp(:);
end

end
