function [S] = EstimateSingularityWithBatches(pointsSet, m, kernel,delta)

n = size(pointsSet, 1);
% bins = false(n,n); % which points I stiil need to check

%% Batch info
nBatched = 4;
skip = floor(n / nBatched);
perm = randperm(n);
bins = false(n,nBatched);

%% Batch Loop
for iBatch = 1 : nBatched
    % run over u to find x in S
    indcesForBatch = perm((iBatch - 1) * skip + 1 :  iBatch * skip);
    tmpPointSetForBatch = pointsSet(indcesForBatch, :);
    
    nPointInBatch = length(indcesForBatch);
    binsBatch = false(n,nPointInBatch); % which points I stiil need to check
      
    %% Check LOOP
    parfor iPoint = 1 : nPointInBatch

        u = tmpPointSetForBatch(iPoint, :);

        % compute std and mean
        weihgts = kernel(pointsSet(:,1), u(1),pointsSet(:,2), u(2), delta);
        meanW = mean(weihgts);
        stdW = std(weihgts);

        % compute MLS
        [~, appFunc] = MLS2D(pointsSet, m, kernel,delta, u);

        epsilon_u = ComputeEpsiolon(pointsSet, appFunc, weihgts, 3 * stdW);

        % find x
        tmplist = 1 : n;
    %     tmplist = tmplist(bins);
        tmplist = tmplist(weihgts(tmplist) > meanW); % first condition
        tmpBin = false(n,1);

        diff = abs(appFunc(pointsSet(tmplist, 1).',pointsSet(tmplist, 2).').' - pointsSet(tmplist, 3));
        tmpBin(diff > (2 * epsilon_u)) = true;

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



        binsBatch(:,iPoint) = tmpBin;
    end
    bins(:, iBatch) = any(binsBatch);
end
bins1 = any(bins); %check if any scenario was fit
S = pointsSet(bins1,:);
end


function [epsilon_u] = ComputeEpsiolon(pointsSet, appFunc, weihgts, stdW)

    indBin = weihgts  > stdW;
    n = size(pointsSet, 1);
    ind = 1 : n;
    ind = ind(indBin);
    
    % Find Max difference
    epsilon_u = max(abs( pointsSet(ind,3) - ...
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