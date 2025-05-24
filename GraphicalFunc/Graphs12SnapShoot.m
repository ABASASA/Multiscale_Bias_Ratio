function [figs] = Graphs12SnapShoot(nDim, testPoints, multiSanpshoot, sourceSanpshoot, singleSanpshoot)
figs = {};
if nDim == 2
    figs = SnapShoot2D(testPoints, multiSanpshoot, sourceSanpshoot, singleSanpshoot);
else
    error('Need to implement');
end
end

function [figs] = SnapShoot2D(testPoints, multiSanpshoot, sourceSanpshoot, singleSanpshoot)
figs = {};
sizeSq = sqrt(length(testPoints));
reshpeFunc = @(A) reshape(A, [sizeSq,sizeSq]);


fig = figure;
subplot(3,1,1);
imagesc((reshpeFunc(sourceSanpshoot).'));
title('Source');
colorbar;

subplot(3,1,2);
imagesc((reshpeFunc(multiSanpshoot).'));
title('Multi');
colorbar;

subplot(3,1,3);
imagesc((reshpeFunc(singleSanpshoot).'));
title('Single');
colorbar;

figs = fig;

end