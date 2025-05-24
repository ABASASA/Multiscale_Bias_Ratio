%% This function comb   ine and presents differents saving via  on measuemrs

% addpath('Graphs');
% hs = [0.2792, 0.2227, 0.1719, 0.1232, 0.0882, 0.0753];
% % factorss = [0.2, 0.14802, 0.10954, 0.081072, 0.075 0.06, 0.05];
% factorss = round(logspace(log10(0.2),log10(0.025),7), 2);
% factorss = factorss(1:6);
% xs =  (1 ./ factorss).^2; % we use it to compute the number of points for each factor
% param = 0.075;
% funcPath = @(numInStr) ['Graphs\280225_Exp_MLS_dist_' num2str(numInStr) '\data.mat'];
% savingPath = ['Graphs\280225_Exp_MLS_Unite\'];
% isSPDFlag = false; % This flag is just for minor grpahicks adjucjemts
%                  
%% Initlize structs
leg = {};
meanMSEsMulti = zeros(length(factorss), 1); %mean -  MSE Bias,var
meanMSEsSingle = zeros(length(factorss), 1); %mean -  MSE Bias,var

upperMSEsMulti = zeros(length(factorss), 1); 
upperMSEsSingle = zeros(length(factorss), 1); 

lowerMSEsMulti = zeros(length(factorss), 1); 
lowerMSEsSingle = zeros(length(factorss), 1); 

meanMSEsMultiRatio = zeros(length(factorss), 1); %mean -  MSE Bias,var
meanMSEsSingleRatio = zeros(length(factorss), 1); %mean -  MSE Bias,var

upperMSEsMultiRatio = zeros(length(factorss), 1); %mean -  MSE Bias,var
upperMSEsSingleRatio = zeros(length(factorss), 1); %mean -  MSE Bias,var

lowerMSEsMultiRatio = zeros(length(factorss), 1); %mean -  MSE Bias,var
lowerMSEsSingleRatio = zeros(length(factorss), 1); %mean -  MSE Bias,var



%% Extract data
for iFactor = 1 : length(factorss)
    path = funcPath(num2str(factorss(iFactor)));
    dataa = load(path);
    tradeoffDeltaMISE1 = dataa.tradeoffDeltaMISE1;
    A = 100 .* tradeoffDeltaMISE1(end, :, 2) ./ ...
            tradeoffDeltaMISE1(end, :, 1);
        
    tradeoffDeltaMISE2 = dataa.tradeoffDeltaMISE2;
    B = 100 .* tradeoffDeltaMISE2(end, :, 2) ./ ...
            tradeoffDeltaMISE2(end, :, 1);
       
     
%      %% MSE - By STD
%      % Mean
%      meanMSEsMulti(iFactor) = mean(tradeoffDeltaMISE1(end, :, 1));
%      meanMSEsSingle(iFactor) = mean(tradeoffDeltaMISE2(end, :, 1));
%      
%      % MSE -  Pre-compute STD
%      tmpStdMulti = std(tradeoffDeltaMISE1(end, :, 1));
%      tmpStdSingle = std(tradeoffDeltaMISE2(end, :, 1));
%      
%      % MSE - Lower bound by STD
%      lowerMSEsMulti(iFactor) = max(meanMSEsMulti(iFactor) - tmpStdMulti, 0); 
%      lowerMSEsSingle(iFactor) = max(meanMSEsSingle(iFactor) - tmpStdSingle, 0); 
%      
%      % MSE - Upper bound by STD
%      upperMSEsMulti(iFactor) = meanMSEsMulti(iFactor) + tmpStdMulti; 
%      upperMSEsSingle(iFactor) = meanMSEsSingle(iFactor) + tmpStdSingle; 
%      
%      %% Ratio - By STD
%      % Mean
%      meanMSEsMultiRatio(iFactor) = mean(A);
%      meanMSEsSingleRatio(iFactor) = mean(B);
%      
%      %Ratio - Pre-Compute STD
%      stdMultiTmp = std(A);
%      stdSingleTmp = std(B);
%      
%      % Ratio - Upper by STD
%      upperMSEsMultiRatio(iFactor) = min(100, meanMSEsMultiRatio(iFactor) + stdMultiTmp);
%      upperMSEsSingleRatio(iFactor) = min(100, meanMSEsSingleRatio(iFactor) + stdSingleTmp); %mean -  MSE Bias,var
% 
%      % Ratio - Lower by STD
%      lowerMSEsMultiRatio(iFactor) = max(0, meanMSEsMultiRatio(iFactor) - stdMultiTmp);
%      lowerMSEsSingleRatio(iFactor) = max(0, meanMSEsSingleRatio(iFactor) - stdSingleTmp); %mean -  MSE Bias,var
%     
     
     %% MSE - By Precntile
     % Mean
     meanMSEsMulti(iFactor) = prctile(tradeoffDeltaMISE1(end, :, 1), 50);
     meanMSEsSingle(iFactor) = prctile(tradeoffDeltaMISE2(end, :, 1), 50);
     
     % MSE - Lower bound by prec
     lowerMSEsMulti(iFactor) = max(prctile(tradeoffDeltaMISE1(end, :, 1), 25), 0); 
     lowerMSEsSingle(iFactor) = max(prctile(tradeoffDeltaMISE2(end, :, 1), 25), 0); 
     
     % MSE - Upper bound by prec
     upperMSEsMulti(iFactor) = prctile(tradeoffDeltaMISE1(end, :, 1), 75); 
     upperMSEsSingle(iFactor) = prctile(tradeoffDeltaMISE2(end, :, 1), 75); 
     
     %% Ratio - By prec
     % Mean
     meanMSEsMultiRatio(iFactor) = prctile(A, 50);
     meanMSEsSingleRatio(iFactor) = prctile(B, 50);
         
     % Ratio - Upper by prec
     upperMSEsMultiRatio(iFactor) = min(100,prctile(A, 75));
     upperMSEsSingleRatio(iFactor) = min(100, prctile(B, 75)); 

     % Ratio - Lower by STD
     lowerMSEsMultiRatio(iFactor) = max(0, prctile(A, 25));
     lowerMSEsSingleRatio(iFactor) = max(0, prctile(B, 25)); %mean -  MSE Bias,var
     
end

maxLimRatio = max(max(upperMSEsMultiRatio), max(upperMSEsSingleRatio));
minLimRatio = min(min(lowerMSEsSingleRatio), min(lowerMSEsMultiRatio));

maxLim = max(max(upperMSEsMulti), max(upperMSEsSingle));
minLim = min(min(lowerMSEsMulti), min(lowerMSEsSingle));

%% Plot

fig1 = figure;
% subplot(2,1,1)
%sts for multi
[dataFlip1] = AreaPlot(xs, lowerMSEsMultiRatio', upperMSEsMultiRatio');
[dataFlip2] = AreaPlot(xs, lowerMSEsSingleRatio', upperMSEsSingleRatio');

fill(dataFlip1.X, dataFlip1.Y, 'b', dataFlip2.X, dataFlip2.Y,'r','FaceAlpha',0.1);
hold on;
plot(xs, upperMSEsMultiRatio, 'b--',...
    xs, lowerMSEsMultiRatio, 'b--',...
    xs, upperMSEsSingleRatio, 'r--',...
    xs, lowerMSEsSingleRatio, 'r--')
plot(xs, meanMSEsMultiRatio, 'b*-', xs,meanMSEsSingleRatio, 'r*-');
xlim([xs(1), xs(end)]);
ylim([minLimRatio, maxLimRatio]);

if isSPDFlag
%     ylim([0, maxLimRatio]);
    legendPos = 'Northeast';
else
%     ylim([0,100]);
    legendPos = 'Southwest';
end

ylabel(' Br = Bias\textsuperscript{2} / MSE \fontsize{8}{0} \selectfont [\%]', 'interpreter','latex')
xlabel('Number of Points in Dataset',  'interpreter','latex');


L(1) = plot(nan, nan, 'b-');
L(2) = plot(nan, nan, 'r-');

legend(L, {'Multiscale', 'Singlescale'}, 'Location', legendPos,...
                        'interpreter','latex', 'fontsize',15);

set(gca,'fontsize',15)
set(gca, 'XScale', 'linear');
% subplot(2,1,2);
fig2 = figure;  
[dataFlip1] = AreaPlot(xs, lowerMSEsMulti', upperMSEsMulti');
[dataFlip2] = AreaPlot(xs, lowerMSEsSingle', upperMSEsSingle');

fill(dataFlip1.X, dataFlip1.Y, 'b', dataFlip2.X, dataFlip2.Y,'r','FaceAlpha',0.1);
hold on;
semilogy(xs, lowerMSEsMulti, 'b--',...
     xs, upperMSEsMulti, 'b--',...
     xs, upperMSEsSingle, 'r--',...
     xs, lowerMSEsSingle, 'r--')

semilogy(xs, meanMSEsMulti, 'b*-', xs, meanMSEsSingle, 'r*-');
ylabel('MSE', 'interpreter','latex')
xlabel('Number of Points in Dataset', 'interpreter','latex')
xlim([xs(1), xs(end)]);

if isSPDFlag
    ylim([minLim, maxLim]);
%     yticks([0.005,0.01,0.02])
    legendPos = 'Southwest';
else
    ylim([minLim, maxLim]);
    yticks([0.00001, 0.0001, 0.001,0.01, 0.1])
    yticklabels({'10^{-5}','10^{-4}','10^{-3}', '10^{-2}', '10^{-1}'})
    legendPos = 'Southwest';
end


L(1) = plot(nan, nan, 'b-');
L(2) = plot(nan, nan, 'r-');
legend(L, {'Multiscale', 'Singlescale'},'Location',legendPos, 'interpreter','latex' ,'fontsize',15)
set(gca,'fontsize',15);
set(gca, 'YScale', 'log');
set(gca, 'XScale', 'linear');
%% Save figures
disp(['Hey, You are running my code. Your Saveing Path is ', savingPath]);
disp('If you want to save the figure Press: 1 else 0');
toSave = input('Enter 1 or 0: ');
if(toSave)
    mkdir(savingPath);
    fileName = ['Br'];
    saveas(fig1,[savingPath, fileName, '.fig']);
    saveas(fig1,[savingPath, fileName, '.jpg']);
    saveas(fig1,[savingPath, fileName, '.pdf']);
    
    
    fileName = ['MSE'];
    saveas(fig2,[savingPath, fileName, '.fig']);
    saveas(fig2,[savingPath, fileName, '.jpg']);
    saveas(fig2,[savingPath, fileName, '.pdf']);
end


%% Compare the two fittings
AMulti = polyfit(log10(xs),log10(meanMSEsMulti'),1);
ASingle = polyfit(log10(xs),log10(meanMSEsSingle'),1);

multiLine = log10(xs) * AMulti(1) + AMulti(2);
singleLine = log10(xs) * ASingle(1) + ASingle(2);

errorMeanMulti = mean(multiLine - log10(meanMSEsMulti'));
R2Multi = 1 - sum((multiLine - log10(meanMSEsMulti')).^2) / sum((multiLine - errorMeanMulti).^2);

errorMeanSingle = mean(singleLine - log10(meanMSEsSingle'));
R2Single = 1 - sum((singleLine - log10(meanMSEsSingle')).^2) / sum((singleLine - errorMeanSingle).^2);

figure;
subplot(2,1,1);
scatter(log10(xs),log10(meanMSEsMulti'));
hold on;
plot(log10(xs),multiLine)
title(['Multiscale: R^2 = ' num2str(R2Multi) ' -- M = ' num2str(AMulti(1))]);
xlabel('log_{10}(Fill distance (h))')
ylabel('log_{10}(MSE)')

subplot(2,1,2);
scatter(log10(xs),log10(meanMSEsSingle'));
hold on;
plot(log10(xs),singleLine)
title(['Singlescale: R^2 = ' num2str(R2Single) ' -- M = ' num2str(ASingle(1))]);

xlabel('log_{10}(Fill distance (h))')
ylabel('log_{10}(MSE)')