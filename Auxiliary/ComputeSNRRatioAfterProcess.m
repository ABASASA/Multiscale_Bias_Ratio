% assumptions: the data is from [0,1]^2.
% It can process only SPD manifold 
function [SNRRatio] = ComputeSNRRatioAfterProcess(func)

[funcNoised] = AddGuassionNoiseSPD(func, 1, 0, 0);
[pointsSetNoised, ~, ~] = SampleRnfunctionEqui(...
                                    funcNoised, [0,1], 0.01, 2);
[funcClean] = AddGuassionNoiseSPD(func, 0, 0, 0);
[pointsSetClean, ~, ~] = SampleRnfunctionEqui(...
                                    funcClean, [0,1], 0.01, 2);

pointsSetClean = pointsSetClean(:,3:end);
pointsSetNoised = pointsSetNoised(:,3:end) - pointsSetClean;

N = size(pointsSetClean,1);
SNRS = zeros(N, 1);
for i = 1 : N
    
   SNRS(i) = norm(pointsSetClean(i,:),2)^2 ./ norm(pointsSetNoised(i,:),2)^2;
end

SNRRatio = mean(SNRS)
median(mean(SNRS))
std(SNRS)
prctile(SNRS, 25)
prctile(SNRS, 75)

end
