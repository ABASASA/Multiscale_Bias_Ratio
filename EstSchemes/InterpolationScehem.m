function [currentEstimation, additionalInput] = InterpolationScehem(XjMod, delta,...
                    kernel, additionalInput)


    
    % St    art by calculating A
    A = KernelMatrix(XjMod, XjMod, kernel, delta, true);
%     figure;
%     imagesc(A)
    % compute alpha
    alpha = linsolve(A, XjMod(:,3));

    % compute s
    currentEstimation = @(x,y) alpha' * kernel(XjMod(:,1), x, XjMod(:,2), y, delta);
    
end