function [funcNoised] = AddGuassionNoiseSPD(func, sigma, pOut, sigmaOut)

    [funcNoised] = @(x) NoisedFuncV1SPD(x, func, sigma, pOut, sigmaOut);
end

function [values, ps] = NoisedFuncV1SPD(x, func, sigma, pOut, sigmaOut)

    fx = func(x);
    ps = false(size(x,1),1);

%     if sigma == 0 
%         values = fx;
%         return;
%     end
    n = sqrt(size(fx,2));
    values = zeros(size(fx,1), size(fx,2));
    for iPoint = 1 : size(fx,1)
        noise = randn(n,n);
        noise = triu(noise) + triu(noise).' - diag(diag(noise));
%         noise = noise + noise.';
        AA = logm(reshape(fx(iPoint,:),n,n));

        tmpSPD = expm(reshape(AA(:)' .*(1 + sigma .* noise(:).'), n,n));
        A = reshape(tmpSPD,n,n);
        if min(eig(tmpSPD)) <=0 | sum(abs(imag(A))) >0 | norm(A - A','fro')>0
            error('Not SPD');
        end
%         eig(tmpSPD)
%         GeodesicDistanceSPD(fx(iPoint,:), tmpSPD(:))
%         [fx(iPoint,:);  tmpSPD(:)']
        values(iPoint,:) =  tmpSPD(:);
    end
end