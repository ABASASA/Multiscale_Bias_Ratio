function [r] = FixedKernelFunction(x1 , x2, y1, y2, delta)

    r = zeros(size(x1));
    dist = sqrt((x1-x2).^2 + (y1-y2).^2);
    r(dist <= delta) = delta^(-2);
end