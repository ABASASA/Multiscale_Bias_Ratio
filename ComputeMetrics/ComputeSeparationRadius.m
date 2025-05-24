% X, Y - are posistion vectors
% Assuming the space is a square
function [separationRadius, indexOfRemovable] = ComputeSeparationRadius(X, Y)


n = length(X);

distance = inf;
cyclicX = circulant(X);
cyclicX = cyclicX(:,2:end);
cyclicY = circulant(Y);
cyclicY = cyclicY(:,2:end);

dists = sqrt((X - cyclicX).^2 + (Y - cyclicY).^2);
[separationRadius, indexOfRemovable] = min(dists(:)); 
[indexOfRemovable, ~] = ind2sub([size(cyclicX,1), size(cyclicX,2)], indexOfRemovable);
% The index of removable is for the thinning algo "GreedyThiningBroneWende"

separationRadius = separationRadius / 2;

end