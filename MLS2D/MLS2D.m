function [polyEst, appFunc, co] = MLS2D(pointsSet, m, kernel,delta, x)

n = size(pointsSet,1);

% Crearte basis to \Pi_m
% basisFunc1 =  polyBasis('canonical',m,2,{'x','y'});
% basisFunc = @(x1,y1) [1, basisFunc1(x1,y1)];

basisFunc = GeneratePolyBasisCanonical(m);

Q = length(basisFunc(x(1), x(2)));

% Create D
D = diag(kernel(pointsSet, x,delta));

% Create P
% P = zeros(n, Q);
% for i = 1 : size(pointsSet,1)
%     P(i,:) = basisFunc(pointsSet(i,1), pointsSet(i,2));  
% end
P = basisFunc(pointsSet(:,1).', pointsSet(:,2).');  

% f
% warning('off','all')
f = pointsSet(:,3);
C = P'*D*P;
co = cond(C);
A = f' * D * P * inv(C);
polyEst = A * (basisFunc(x(1), x(2))'); % p(0)
appFunc = @(x,y) A * (basisFunc(x, y)'); % (p(x)
% warning('on','all')
end
