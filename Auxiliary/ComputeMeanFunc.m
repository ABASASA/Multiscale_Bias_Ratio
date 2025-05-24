function [mean1] = ComputeMeanFunc(func, gridLimits, dist)

[A, ~] = Sample2DfunctionEqui(func, gridLimits, dist);
mean1 = mean(A(:,3));

end