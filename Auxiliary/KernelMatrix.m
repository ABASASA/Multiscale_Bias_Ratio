
function [r] = KernelMatrix(Xj, Xk, kernel, delta, flag)
[X1, X2] = meshgrid(Xj(:,1), Xk(:,1));
[Y1, Y2] = meshgrid(Xj(:,2), Xk(:,2));
r = kernel(X1, X2, Y1, Y2, delta);
% figure
% imagesc(r)
if flag
    disp(['Condition of kerenl: ', num2str(cond(r))]);
end
% Test Via element wise
% r1 = zeros(size(r));

% for i = 1: size(Xj,1)
%     for j = 1: size(Xk,1)
%        r1(i,j) = kernel(Xj(i,1), Xk(j,1), Xj(i,2), Xk(j,2), delta);
%     end
% end
% norm(r1(:) - r(:))
end
