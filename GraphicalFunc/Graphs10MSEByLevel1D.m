function [figs] = Graphs10MSEByLevel1D(scale, lebelsScale, MSE_Multiscale, MSE_Singlescale,...
                         testPoints, nLevels, nDim, tradeoffInegrated1, tradeoffInegrated2...
                         , meanMulti, meanSingle)

 figs = {};
 if nDim == 1
     figs = Graph10_1D(MSE_Multiscale, MSE_Singlescale,...
                         testPoints, nLevels);
 elseif nDim == 2
%      figs = Graph10_1D(MSE_Multiscale, MSE_Singlescale,...
%                          1 : length(testPoints), nLevels, [3]);
     figs = Graph10_2D(scale, lebelsScale, MSE_Multiscale, MSE_Singlescale,...
                         1 : length(testPoints), nLevels, tradeoffInegrated1, tradeoffInegrated2...
                         , meanMulti, meanSingle);
 end
                     
end


function [figs] = Graph10_1D(MSE_Multiscale, MSE_Singlescale,...
                         testPoints, nLevels, indLevel)
if ~exist('indLevel','var')
    indLevel = 1  : nLevels;
end
    
  
 
 figs = {};
 fig = figure;
 colors = zeros(length(indLevel), 3);
 for i = 1 : 3 % 1 - bias 2 - variance, 3 MSE
    subplot(3, 1, i);
    leg = {};

    for iLevel = 1 : length(indLevel)
        
        hold on;
        pp = semilogy(testPoints, MSE_Multiscale(indLevel(iLevel), :, i));
        leg{iLevel} = num2str(indLevel(iLevel));
        colors(iLevel, :) = pp.Color;

    end
   
    for iLevel = 1 : length(indLevel)
        hold on;
        semilogy(testPoints, MSE_Singlescale(indLevel(iLevel), :, i),...
                            '--', 'Color', colors(iLevel, :))
    end
    if i == 1
        legend(leg);
    end
    
    switch i
        case 2
            title('Bias');
        case 3
            title('Variance');
        case 1
            title('MSE');
    end
 end
figs{1} = fig;


%% Just Plot MSE
fig = figure;
leg = {};
 colors = zeros(length(indLevel), 3);

for iLevel = 1 : length(indLevel)
    hold on;
    pp = semilogy(testPoints, MSE_Multiscale(indLevel(iLevel), :, i));
    leg{iLevel} = num2str(indLevel(iLevel));
    colors(iLevel, :) = pp.Color;
    hold on;
end
semilogy(testPoints, MSE_Singlescale(end, :, 1), '--', 'Color', colors(end, :))
leg{iLevel + 1} = 'Single';
legend(leg);
title('MSE');

figs{2} = fig;

%% Show the change in bias variance tradeoff
fig = figure;
for iLevel = 1 : length(indLevel)
    hold on;
    pp = semilogy(testPoints, 100 .* MSE_Multiscale(indLevel(iLevel), :, 2) ./ ...
            MSE_Multiscale(indLevel(iLevel), :, 1), 'Color', colors(iLevel, :));
    leg{iLevel} = num2str(indLevel(iLevel));
    
end
semilogy(testPoints, 100 .* MSE_Singlescale(end, :, 2) ./ ...
            MSE_Singlescale(end, :, 1), '--', 'Color', colors(end, :));
leg{iLevel + 1} = 'Single';

ylabel('Bias From MSE [%]');
legend(leg);
figs{3} = fig;
end


function figs = Graph10_2D(scale, lebelsScale, MSE_Multiscale, MSE_Singlescale,...
                         testPoints, nLevels, tradeoffInegrated1, tradeoffInegrated2,...
                         meanMulti, meanSingle, indLevel)

% If one wants to present the results for a specific levels use 'indLevel =
% [1,2,4] for example
if ~exist('indLevel','var')
    indLevel = 1  : nLevels;
end
nInRow = sqrt(size(testPoints,2));
figs = {};


%% LEP of data
% Initizlize
ratios = zeros(length(indLevel), size(MSE_Multiscale,2));

meanMSEsMulti = zeros(length(indLevel), 1); 
meanMSEsSingle = zeros(length(indLevel), 1); 

upperMSEsMulti = zeros(length(indLevel), 1); 
upperMSEsSingle = zeros(length(indLevel), 1); 

lowerMSEsMulti = zeros(length(indLevel), 1); 
lowerMSEsSingle = zeros(length(indLevel), 1); 

meanMSEsMultiRatio = zeros(length(indLevel), 1);
meanMSEsSingleRatio = zeros(length(indLevel), 1); 

upperMSEsMultiRatio = zeros(length(indLevel), 1); 
upperMSEsSingleRatio = zeros(length(indLevel), 1);

lowerMSEsMultiRatio = zeros(length(indLevel), 1);
lowerMSEsSingleRatio = zeros(length(indLevel), 1); 

% Extract data
for iLevel = 1 : length(indLevel)
    A =  100 .* MSE_Multiscale(indLevel(iLevel), :, 2) ./ ...
            MSE_Multiscale(indLevel(iLevel), :, 1);
    B =  100 .* MSE_Singlescale(indLevel(iLevel), :, 2) ./ ...
            MSE_Singlescale(indLevel(iLevel), :, 1);
    ratios(iLevel, :) = A ./ B;
    
    % STD Approach
    meansMulti(iLevel, 1) = mean( A);
    meansMulti(iLevel, 2) = mean( MSE_Multiscale(indLevel(iLevel), :, 1));

    meansSingle(iLevel, 1) = mean( B);
    meansSingle(iLevel, 2) = mean( MSE_Singlescale(indLevel(iLevel), :, 1));
    
    
    stdsMulti(iLevel, 1) = std(A);
    stdsMulti(iLevel, 2) = std( MSE_Multiscale(indLevel(iLevel), :, 1));
    
    stdsSingle(iLevel, 1) = std(B);
    stdsSingle(iLevel, 2) = std( MSE_Singlescale(indLevel(iLevel), :, 1));

      
%      %% MSE - By STD
%      % Mean
%      meanMSEsMulti(iLevel) = mean(MSE_Multiscale(indLevel(iLevel), :, 1));
%      meanMSEsSingle(iLevel) = mean(MSE_Singlescale(indLevel(iLevel), :, 1));
%      
%      % MSE -  Pre-compute STD
%      tmpStdMulti = std(MSE_Multiscale(indLevel(iLevel), :, 1));
%      tmpStdSingle = std(MSE_Singlescale(indLevel(iLevel), :, 1));
%      
%      % MSE - Lower bound by STD
%      lowerMSEsMulti(iLevel) = max(meanMSEsMulti(iLevel) - tmpStdMulti, 0); 
%      lowerMSEsSingle(iLevel) = max(meanMSEsSingle(iLevel) - tmpStdSingle, 0); 
%      
%      % MSE - Upper bound by STD
%      upperMSEsMulti(iLevel) = meanMSEsMulti(iLevel) + tmpStdMulti; 
%      upperMSEsSingle(iLevel) = meanMSEsSingle(iLevel) + tmpStdSingle; 
%      
%      %% Ratio - By STD
%      % Mean
%      meanMSEsMultiRatio(iLevel) = mean(A);
%      meanMSEsSingleRatio(iLevel) = mean(B);
%      
%      %Ratio - Pre-Compute STD
%      stdMultiTmp = std(A);
%      stdSingleTmp = std(B);
%      
%      % Ratio - Upper by STD
%      upperMSEsMultiRatio(iLevel) = min(100, meanMSEsMultiRatio(iLevel) + stdMultiTmp);
%      upperMSEsSingleRatio(iLevel) = min(100, meanMSEsSingleRatio(iLevel) + stdSingleTmp); %mean -  MSE Bias,var
% 
%      % Ratio - Lower by STD
%      lowerMSEsMultiRatio(iLevel) = max(0, meanMSEsMultiRatio(iLevel) - stdMultiTmp);
%      lowerMSEsSingleRatio(iLevel) = max(0, meanMSEsSingleRatio(iLevel) - stdSingleTmp); %mean -  MSE Bias,var
%     
     
     %% MSE - By Precntile
     % Mean
     meanMSEsMulti(iLevel) = prctile(MSE_Multiscale(indLevel(iLevel), :, 1), 50);
     meanMSEsSingle(iLevel) = prctile(MSE_Singlescale(indLevel(iLevel), :, 1), 50);
     
     % MSE - Lower bound by prec
     lowerMSEsMulti(iLevel) = max(prctile(MSE_Multiscale(indLevel(iLevel), :, 1), 25), 0); 
     lowerMSEsSingle(iLevel) = max(prctile(MSE_Singlescale(indLevel(iLevel), :, 1), 25), 0); 
     
     % MSE - Upper bound by prec
     upperMSEsMulti(iLevel) = prctile(MSE_Multiscale(indLevel(iLevel), :, 1), 75); 
     upperMSEsSingle(iLevel) = prctile(MSE_Singlescale(indLevel(iLevel), :, 1), 75); 
     
     %% Ratio - By prec
     % Mean
     meanMSEsMultiRatio(iLevel) = prctile(A, 50);
     meanMSEsSingleRatio(iLevel) = prctile(B, 50);
         
     % Ratio - Upper by prec
     upperMSEsMultiRatio(iLevel) = min(100,prctile(A, 75));
     upperMSEsSingleRatio(iLevel) = min(100, prctile(B, 75)); 

     % Ratio - Lower by STD
     lowerMSEsMultiRatio(iLevel) = max(0, prctile(A, 25));
     lowerMSEsSingleRatio(iLevel) = max(0, prctile(B, 25)); 
      
    
end
% Limits
maxLimRatio = max(max(upperMSEsMultiRatio), max(upperMSEsSingleRatio));
minLimRatio = min(min(lowerMSEsSingleRatio), min(lowerMSEsMultiRatio));

maxLim = max(max(upperMSEsMulti), max(upperMSEsSingle));
minLim = min(min(lowerMSEsMulti), min(lowerMSEsSingle));

% First figure RAtio
fig = figure;

[dataFlip1] = AreaPlot(scale, lowerMSEsMultiRatio', upperMSEsMultiRatio');
[dataFlip2] = AreaPlot(scale, lowerMSEsSingleRatio', upperMSEsSingleRatio');

fill(dataFlip1.X, dataFlip1.Y, 'b', dataFlip2.X, dataFlip2.Y,'r','FaceAlpha',0.1);
hold on;
plot(scale, upperMSEsMultiRatio, 'b--',...
    scale, lowerMSEsMultiRatio, 'b--',...
    scale, upperMSEsSingleRatio, 'r--',...
    scale, lowerMSEsSingleRatio, 'r--')
plot(scale, meanMSEsMultiRatio, 'b*-', scale,meanMSEsSingleRatio, 'r*-');
xlim([scale(1), scale(end)]);
ylim([minLimRatio, maxLimRatio]);
% ylim([0,100])
% ylim([0, maxLimRatio]);

ylabel(' Br = Bias\textsuperscript{2} / MSE \fontsize{8}{0} \selectfont [\%]', 'interpreter','latex')
xlabel('Data [\%]', 'interpreter','latex');
% ylim([0,100]);
xticks([70,80,90,100]);
L(1) = plot(nan, nan, 'b-');
L(2) = plot(nan, nan, 'r-');
legend(L, {'Multiscale', 'Singlescale'}, 'Location','southwest', 'interpreter','latex', 'fontsize',15)
set(gca,'fontsize',15);

figs{end + 1} = fig;

% MSE
fig = figure;

[dataFlip1] = AreaPlot(scale, lowerMSEsMulti', upperMSEsMulti');
[dataFlip2] = AreaPlot(scale, lowerMSEsSingle', upperMSEsSingle');

fill(dataFlip1.X, dataFlip1.Y, 'b', dataFlip2.X, dataFlip2.Y,'r','FaceAlpha',0.1);
hold on;
semilogy(scale, lowerMSEsMulti, 'b--',...
     scale, upperMSEsMulti, 'b--',...
     scale, upperMSEsSingle, 'r--',...
     scale, lowerMSEsSingle, 'r--')

semilogy(scale, meanMSEsMulti, 'b*-', scale, meanMSEsSingle, 'r*-');
ylabel('MSE', 'interpreter','latex')
xlabel('Data [\%]', 'interpreter','latex')
xlim([scale(1), scale(end)]);
ylim([minLim, maxLim]);
% ylim([0.001, maxLim]);
% yticks([0.005,0.01,0.02])
% yticklabels({'5\cdot10^{-3}', '10^{-2}', '2\cdot10^{-2}'})
xticks([70,80,90,100]);
L(1) = plot(nan, nan, 'b-');
L(2) = plot(nan, nan, 'r-');
legend(L, {'Multiscale', 'Singlescale'},'Location','northeast', 'interpreter','latex' ,'fontsize',15)
set(gca,'fontsize',15);
% set(gca, 'YScale', 'log');

figs{end + 1} = fig;

%%
% fig = figure;
% leg = {};
% for iLevel = 1: length(indLevel)
%     hold on;
%     plot(cet, LEPsSinge(indLevel(iLevel),:) ./ LEPS(indLevel(iLevel),:), '*-');
%     leg{end + 1} = num2str(indLevel(iLevel));
% end
% hold on;
% title('Ratios - Bigger than 1 is goog')
% legend(leg)
% % plot([cet(1), cet(end)], [1,1], 'k')
% xlabel('Precentile');
% title('Ratio of per(Br(Single),p) / per(Br(Multi),p)')
% figs{end + 1} = fig;

%% Box Plot
% fig = figure;
% ratios = 1./ratios;
% fig = BoxPlotAsaf(fig, scale, ratios, 'b*-');
% ylim([0, 5]);
% ylabel('Ratio of Br - Single / Multiscale');
% xlabel('Precentages of the data');
% 
% title('Ratios Br(Single) / Br(Multiscale) ')
% figs{end + 1} = fig;
%% MISE, IBias and  IVar - MISE(right plot)
fig = figure;
ratioMultiscale = 100 * (tradeoffInegrated1(indLevel,2) ./ tradeoffInegrated1(indLevel,1));
ratioSinglescale = 100 * (tradeoffInegrated2(indLevel,2) ./ tradeoffInegrated2(indLevel,1));
subplot(2,1,1);
yyaxis left
plot(scale, ratioMultiscale, '*-', scale, ratioSinglescale, '*--');
xlabel(lebelsScale)
ylabel('IBias^2 / MISE [%]')

% subplot(2,1,2);
yyaxis right

plot(scale, log10(tradeoffInegrated1(:,1)), '*-', scale, log10(tradeoffInegrated1(:,2)), '*--');
ylabel('log_{10}(MISE)')
xlabel(lebelsScale)
figs{end + 1} = fig;

legend('Multiscale', 'Singlescale');

%% MISE, IBias and  IVar - mean (right plot)

subplot(2,1,2);
yyaxis left
plot(scale, ratioMultiscale, '*-', scale, ratioSinglescale, '*--');
xlabel(lebelsScale)
ylabel('IBias^2 / MISE [%]')

% subplot(2,1,2);
yyaxis right
plot(scale, log10(meanMulti), '*-', scale, log10(meanSingle), '*--');
ylabel('log_{10}(mean error)')
xlabel(lebelsScale)
figs{end + 1} = fig;

legend('Multiscale', 'Singlescale');



%%
% 
% reShFunc = @(A) reshape(A,[nInRow, nInRow]);
% for iLevel = 1 : length(indLevel)
%     fig = figure;
%     hold on;
%     subplot(2,1,1);
%     A = reShFunc(100 .* MSE_Multiscale(indLevel(iLevel), :, 2) ./ ...
%             MSE_Multiscale(indLevel(iLevel), :, 1));
%     imagesc(A.');
%     title(['Level ' num2str(indLevel(iLevel)) ': Multscale - Bias From MSE [%]']);
% 
%     colorbar;
%     caxis([0 100]);
% 
%     subplot(2,1,2);
%     B = reShFunc(100 .* MSE_Singlescale(indLevel(iLevel), :, 2) ./ ...
%             MSE_Singlescale(indLevel(iLevel), :, 1));
%     imagesc(B.');
%     colorbar;
%     caxis([0 100]);
%     title(['Level ' num2str(indLevel(iLevel)) ': Singlescale - Bias From MSE [%]']);
% 
%     figs{end + 1} = fig;
% end


end
