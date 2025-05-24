function [fig, max1, max2] = Graphs1DEnvelopeNew(levels, xlab, errorsDelta1, errorsDelta2)

nLevels = size(errorsDelta1,2);

mean1 = squeeze(mean(errorsDelta1(:,:, 1), 1));

median1 = squeeze(mean(errorsDelta1(:,:, 2), 1));

max1 = squeeze(mean(errorsDelta1(:,:, 3), 1));


mean2 = squeeze(mean(errorsDelta2(:,:, 1), 1));

median2 = squeeze(mean(errorsDelta2(:,:, 2), 1));

max2 = squeeze(mean(errorsDelta2(:,:, 3), 1));

%% Plot error
ylab = 'log_{10}(error)';
leg = {'Multiscale', 'Single'};
% xlab = 'Scale';
% levels = scale;
fig = figure;

% one = ones(nLevels,1);
subplot(3,1,1)
hold on;
plot(levels, log10(mean1),'-*', levels, log10(mean2) , '-*');
title('Mean Mean error')
xlabel(xlab)
ylabel(ylab)
legend(leg);

subplot(3,1,2)
hold on;
plot(levels, log10(median1),'-*',levels, log10(median2),'-*');
title('Mean Median error')
ylabel(ylab)
xlabel(xlab)

subplot(3,1,3)
hold on;
plot(levels, log10(max1),'-*',levels, log10(max2), '-*');
title('Mean Max error')
ylabel(ylab)
xlabel(xlab)

end