function [S] = EstimateSingularity(pointsSet, m, kernel,delta)

n = size(pointsSet, 1);
bins = false(n,n); % which points I stiil need to check
% run over u to find x in S
parfor iPoint = 1 : n
    
    u = pointsSet(iPoint, :);
    
    % compute std and mean
    weihgts = kernel(pointsSet(:,1), u(1), pointsSet(:,2), u(2), delta);
    meanW = mean(weihgts);
    stdW = std(weihgts);
    
    %% compute MLS
    [~, appFunc] = MLS2D(pointsSet, m, kernel,delta, u);
    
    %% Compute epsilon_u
    epsilon_u = ComputeEpsiolon(pointsSet, appFunc, weihgts, 3 * stdW);
    
    %% find which x are fit for this u
    tmplist = 1 : n;
    tmplist = tmplist(weihgts > meanW); % first condition
    tmpBin = false(n,1);
    
    diff = -inf(n,1);
    diff(tmplist) = abs(appFunc(pointsSet(tmplist, 1).',pointsSet(tmplist, 2).').'...
                            - pointsSet(tmplist, 3));
    tmpBin(diff > (2 * epsilon_u)) = true;
    
    bins(:,iPoint) = tmpBin;
%     if (sum(tmpBin) > 0)
%         disp(num2str(iPoint));
%     end
%     tmplist = 755;
%     tmplist = tmplist(weihgts(tmplist) > meanW); % first condition
%     diff = abs(appFunc(pointsSet(tmplist, 1).',pointsSet(tmplist, 2).').'...
%                             - pointsSet(tmplist, 3));
%     A = diff >= (2 * epsilon_u);
%     if A
%         disp(755);
%     end
    % OLD code in loop
%     tmpBin1 = false(n,1);
% 
%     for j = 1 : length(tmplist)
%         tmpPoint = pointsSet(tmplist(j),:);
%         diff = abs(appFunc(tmpPoint(1),tmpPoint(2)) - tmpPoint(3));
%         if ( diff > (2 * epsilon_u))
%            tmpBin(tmplist(j)) = true;
%         end
%     end
    
end

bins1 = (sum(bins,2) ~= 0); %check if any scenario was fit
S = pointsSet(bins1,:);
end


function [epsilon_u] = ComputeEpsiolon(pointsSet, appFunc, weihgts, stdW)

    n = size(pointsSet, 1);
    ind = 1 : n;
    ind = ind(weihgts  > stdW);
    
    % Find Max difference
    epsilon_u = max(abs(pointsSet(ind,3) - ...
                appFunc(pointsSet(ind,1).',pointsSet(ind,2).').'));

    % old Code
%     tmpMax = 0;
%     for i = 1 : length(ind)
%         tmpPoint = pointsSet(ind(i),:);
%         
%         tmpMax = max(tmpMax, abs(tmpPoint(3) - appFunc(tmpPoint(1),tmpPoint(2))));
%     end
%     epsilon_u = tmpMax;
    
    
end