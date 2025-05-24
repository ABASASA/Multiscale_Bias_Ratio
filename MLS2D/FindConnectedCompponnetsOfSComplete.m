function [V, comp] = FindConnectedCompponnetsOfSComplete(pointsSet, S, tresh1, tresh2)

N = size(pointsSet, 1);
numS = size(S, 1);

%% start by finding the vertices in the graph 
VFlags = false(N, 1);

for i  = 1 : numS
    tmpPoint = S(i, :);
    dists = sqrt((pointsSet(:,1) - tmpPoint(1)).^2 + (pointsSet(:,2) - tmpPoint(2)).^2);
    
    VFlags(dists <= tresh1) = true;
    
    if sum(VFlags) + numS == N
        break;
    end
    
end

V = pointsSet(~VFlags, :); % the set of points such that it is not close to the hyperplane
 


%% Find the edges
% pdist compute all the sitances between the points.
% squareform - turn it to matrix (0 in the diagonal)
% The condition (<= tresh2) test if the distance is small enough also deal
% with the 0 in the digonal
E = squareform(pdist(V(:,1:2))) <= tresh2;

%% compute components

g = graph(E);
comp  = conncomp(g);


end