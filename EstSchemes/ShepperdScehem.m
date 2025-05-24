function [currentEstimation, additionalInput] = ShepperdScehem(XjMod, delta,...
                    kernel, additionalInput, nDim)

    currentEstimation = @(X) ShepperdFunctionV1(X, XjMod , delta,...
                    kernel, nDim);
    
end

function [est] = ShepperdFunctionV1(X, XjMod, delta,...
                    kernel, nDim)
    
    tmp = kernel(X, XjMod(:,1:nDim), delta);
    a = 1 - tmp;
    a(a<0) = 0;
    tmp = a.^4 .* (4 .* tmp + 1);
    est = XjMod(:,nDim + 1)' * tmp / sum(tmp);
                        
end

% Here delta = r
% According to Shepperd paper from 1968
function [est] = ShepperdFunctionV2(x,y, XjMod, delta,...
                    kernel)
    ind1 = x == XjMod(:,1);
    ind2 = y == XjMod(:,2);
    if sum(ind1 & ind2) > 0 
       est = XjMod(ind1 & ind2, 3);
       return;
    end
           
    dists = sqrt((x-XjMod(:,1)).^2 + (y-XjMod(:,2)).^2);
    
    nCp = sum(dists <= delta);
    rp = 0; % the r tag for P
    if nCp < 5
        distsS = sort(dists);
        rp = distsS(5);
    elseif nCp < 11
        rp = delta;
    else
        distsS = sort(dists);
        rp = distsS(11);
    end
    
    weights = zeros(length(dists),1);
    tmpDists = dists(dists <= rp/3);
    weights(dists <= rp/3) = 1 ./ tmpDists;
    tmpDists = dists(dists > rp/3 & dists <= rp);
    weights(dists > rp/3 & dists <= rp) = (27 ./ (4 .* rp)) .* (tmpDists ./ rp - 1 ).^2;
    weights(dists > rp) = 0;
    
    est = XjMod(:,3)' * (weights.^2) / sum((weights.^2));
                        
end


