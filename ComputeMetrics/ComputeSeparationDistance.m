% X, Y - are posistion vectors
function [distance] = ComputeSeparationDistance(X, Y)

n = length(X);

Xs = zeros(n, n-1);
Ys = zeros(n, n-1);

for i = 1 : n-1
   
    Xs(:, i) = circshift(X, i);
    Ys(:, i) = circshift(Y, i);
    
end
allDist =  sqrt((X - Xs).^2 + (Y - Ys).^2);

distance = min(allDist(:));

end