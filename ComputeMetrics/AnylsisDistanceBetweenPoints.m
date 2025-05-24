function [] = AnylsisDistanceBetweenPoints(pointSet)

X = pointSet(:,1);
Y = pointSet(:,2);

ndataPoints = 15;
n = size(pointSet,1);

distance = inf;
cyclicX = circulant(X);
cyclicX = cyclicX(:,2:end);
cyclicY = circulant(Y);
cyclicY = cyclicY(:,2:end);

dists = sqrt((X - cyclicX).^2 + (Y - cyclicY).^2);
distsVec = dists(:);

rs = logspace(log10(min(distsVec)), log10(max(distsVec)),ndataPoints);
% compute the number of distances which is smaller than rs(i)
nums = zeros(ndataPoints, 3);
for i = 1 : ndataPoints
   tmp = sum(dists <= rs(i),2);
   nums(i,1) =  mean(tmp);
   nums(i,2) =  max(tmp);
   nums(i,3) =  min(tmp);
   
end
%% plot
figure;
loglog((rs), nums(:,1), '*--',(rs), nums(:,2), '*--',(rs), nums(:,3), '*--')
legend('mean','max','min')
xlabel('log_{10} radii');
ylabel('log_{10} number of neighbers in this distance');

end