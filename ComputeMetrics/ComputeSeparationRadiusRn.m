% X, Y - are posistion vectors
% Assuming the space is a square
function [separationRadius, indexOfRemovable] =...
    ComputeSeparationRadiusRn(samplingPlaces)


n = size(samplingPlaces, 1);
nDim = size(samplingPlaces, 2);

summ = 0;
for j = 1 : nDim
    cyclicJ = circulant(samplingPlaces(:, j));
    cyclicJ = cyclicJ(:,2:end); % Removing the first line (which is equal to hte comparsasion)
    summ = summ + (samplingPlaces(:, j) - cyclicJ) .^2;
end

dists = sqrt(summ);

[separationRadius, indexOfRemovable] = min(dists(:)); 
[indexOfRemovable, ~] = ind2sub([size(cyclicJ,1), size(cyclicJ,2)], indexOfRemovable);
% The index of removable is for the thinning algo "GreedyThiningBroneWende"

separationRadius = separationRadius / 2;

end