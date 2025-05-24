function [compOfX, flag] = ExpansionOfComponnets(V, comp, kernel,...
                                    m, deltas, tresh1, x, fx)
Ncomp = length(unique(comp));
isPossibleInComp = false(Ncomp,1);
values = inf(Ncomp,1);
for i = 1 : Ncomp
        
    tmpSet = V(comp==i,:);
    % check if it is in the conponent
    if sum((sum(abs(tmpSet(:,1:2)-x),2) ==0 )) ~= 0
        compOfX = i;
        flag = true;
        return;
    end
    %% compute the estimation error from the value
    [~, appFunc, co] = MLS2D(tmpSet, m, kernel,deltas, x);
    if co < tresh1
        isPossibleInComp(i) = true;
        values(i) = appFunc(x(1), x(2)) - fx;
    end
end

[~,i] = min(values);

if isPossibleInComp(i)
    compOfX = i;
    flag = true;
else
    flag = false;
    compOfX = [0];
%     error('There is a bug here');
end

end

