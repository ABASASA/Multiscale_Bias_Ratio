function [r] = NormWeighForMLS2DExp(X1, X2, delta, nDim)

%     r = delta^(-1) * (sqrt((x1-x2).^2 + (y1-y2).^2));
    r = 0;
    for i = 1 : nDim
        r = r + (X1(:,i) - X2(:, i)).^2;
    end
    
%     r =  exp(-2 .* (delta^(-1) .* r).^2);
    r =  exp(-2 .* delta^(-2) .* r);

end