function [currentEstimation, additionalInput] = ShepperdSPDScehem(XjMod, delta,...
                    kernel, additionalInput, nDim)

    currentEstimation = @(X) ShepperdFunctionV1(X, XjMod , delta,...
                    kernel, nDim);
    
end

function [est] = ShepperdFunctionV1(X, XjMod, delta,...
                    kernel, nDim)
    nDimOutPut = size(XjMod,2) - nDim;
    %% First, we compute the alphas coeffients
    dists = kernel(X, XjMod(:,1:nDim), delta);
    a = 1 - dists;
    a(a<0) = 0;
    alphas = a.^4 .* (4 .* dists + 1); % Distance with raidl function
    rellevantInd = alphas ~= 0; % Those are the indeces of the rellevant data points
    alphasRellevant = alphas(rellevantInd) ./ sum(alphas);
    alphasRellevantMat = repmat(alphasRellevant, [1, nDimOutPut]);
    
    numOfRelOb = length(alphasRellevant);

    indecesOb = 1:size(XjMod,1);
    indecesOb = indecesOb(rellevantInd);
    
    
    % A small unit test
    if size(alphasRellevantMat, 1) ~= numOfRelOb |...
            size(alphasRellevantMat, 2) ~= nDimOutPut
        error('Problem Dim');
    end

    
    %% initlize the iteration
    Niterations = 3; % a pre-deterrmaint number of iterations (in the fuature it can change to a stopting condition);
    thetas = zeros(Niterations + 1, 1); % The weights of each iterative step
    Ys = zeros(Niterations + 1, nDimOutPut);
    

    Yold = sum(alphasRellevantMat .* XjMod(indecesOb, nDim + 1 : end), 1); 
    Ys(1, :) = Yold;
    condsOld = ComputeCond(Yold, XjMod(indecesOb, :), nDim);
    thetas(1) = ComputeTheta(condsOld, alphasRellevant);
    
    %% Iterative Process
    try 
        for iIter = 1 : Niterations
           logs = LogSPD(XjMod(indecesOb, nDim + 1 : end), Yold);
           tmp =  thetas(iIter) .* sum(alphasRellevantMat .* logs, 1);
           Yold = ExpSPD(tmp, Yold); % compute Guess
           Ys(iIter + 1,  :) = Yold;
    
           if iIter ~= Niterations
             condsOld = ComputeCond(Yold, XjMod(indecesOb, :), nDim);
             thetas(iIter + 1) = ComputeTheta(condsOld, alphasRellevant);
           end
        end
        est = Yold;
    catch
        est =  Ys(1, :);
    end
end

function [conds] = ComputeCond(Ycurrent, XjModReduced, nDim)
%% initlize
conds = zeros(size(XjModReduced, 1), 1); 

numEl = size(XjModReduced, 2) - nDim;
sqNumel = sqrt(numEl);

funcResh = @(X) reshape(X, [sqNumel, sqNumel]); % respahing function

invsqY = inv(sqrtm(funcResh(Ycurrent)));
% compute condition number
for i = 1 : size(XjModReduced, 1)
    conds(i) = cond(invsqY * funcResh(XjModReduced(i, nDim + 1 : end))...
                    * invsqY);
end
end

function [theta] = ComputeTheta(conds, alphas)

if length(conds) ~= length(alphas)
   error('ERror in diemnsions'); 
end

conds = (conds + 1) ./ (conds - 1);
theta = 2 ./ sum(alphas .* conds);

end