%% This function comb   ine and presents differents saving via  on measuemrs
% 
% addpath('Graphs');
% 
% % % R^2
% factorss = [0:0.25:2.5];
% factorss = factorss(2:end);
% 
% param = 0.075;
% funcPath = @(numInStr) ['Graphs\280225_COS_Shp_sigma_'...
%                     num2str(numInStr) '_dist_' num2str(param) '\data.mat'];
% savingPath = ['Graphs\280225_COS_Shp_dist_' num2str(param) '_Unite\'];
% isSPDFlag = false; % This flag is just for minor grpahicks adjucjemts
% 
% % SPD
% factorss = [0.05:0.05:0.25];
% param = 0.075;
% funcPath = @(numInStr) ['Graphs\080325_SPD__newFunc_Shp_sigma_'...
%                     num2str(numInStr) '_dist_' num2str(param) '\data.mat'];
% savingPath = ['Graphs\080325_SPD__newFunc_Shp_Unite_dist_' ...
%                          num2str(param) '_Unite\'];
% isSPDFlag = true; % This flag is just for minor grpahicks adjucjemts
%                  
%% Initlize structs
leg = {};
nParams = length(factorss);
meanMSEsMulti = zeros(nParams, 1); %mean -  MSE Bias,var
meanMSEsSingle = zeros(nParams, 1); %mean -  MSE Bias,var

upperMSEsMulti = zeros(nParams, 1); 
upperMSEsSingle = zeros(nParams, 1); 

lowerMSEsMulti = zeros(nParams, 1);  
lowerMSEsSingle = zeros(nParams, 1); 

meanMSEsMultiRatio = zeros(nParams, 1);  %mean -  MSE Bias,var
meanMSEsSingleRatio = zeros(nParams, 1);  %mean -  MSE Bias,var

upperMSEsMultiRatio = zeros(nParams, 1);  %mean -  MSE Bias,var
upperMSEsSingleRatio = zeros(nParams, 1);  %mean -  MSE Bias,var

lowerMSEsMultiRatio = zeros(nParams, 1);  %mean -  MSE Bias,var
lowerMSEsSingleRatio = zeros(nParams, 1);  %mean -  MSE Bias,var



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
if isSPDFlag
    chagneSNR = factorss;
    for i = 1 : length(factorss)
        chagneSNR(i) = ComputeSNRNumericllySPD(factorss(i));
    end
else
    chagneSNR = 1 ./ factorss.^2;
end
[dataFlip1] = AreaPlot(chagneSNR, lowerMSEsMultiRatio', upperMSEsMultiRatio');
[dataFlip2] = AreaPlot(chagneSNR, lowerMSEsSingleRatio', upperMSEsSingleRatio');

fill(dataFlip1.X, dataFlip1.Y, 'b', dataFlip2.X, dataFlip2.Y,'r','FaceAlpha',0.1);
hold on;
loglog(chagneSNR, upperMSEsMultiRatio, 'b--',...
    chagneSNR, lowerMSEsMultiRatio, 'b--',...
    chagneSNR, upperMSEsSingleRatio, 'r--',...
    chagneSNR, lowerMSEsSingleRatio, 'r--')
loglog(chagneSNR, meanMSEsMultiRatio, 'b*-', chagneSNR, meanMSEsSingleRatio, 'r*-');
% xlim([factorss(1), factorss(end)]);
xlim([chagneSNR(end), chagneSNR(1)]);
ylim([minLimRatio, maxLimRatio]);


if isSPDFlag
    ylim([0, maxLimRatio]);
    legendPos = 'Southeast';
    xticks([2, 5, 10, 40])
    xticklabels({'2', '5', '10', '40'})
else
    ylim([0,100]);
    legendPos = 'Southeast';
    xticks([0.2,1,10])
    xticklabels({'0.2', '1', '10'})
end

ylabel(' Br = Bias\textsuperscript{2} / MSE \fontsize{8}{0} \selectfont [\%]', 'interpreter','latex')
% xlabel('p - Noise to Signal Ratio',  'interpreter','latex');
xlabel('SNR',  'interpreter','latex');


L(1) = plot(nan, nan, 'b-');
L(2) = plot(nan, nan, 'r-');

legend(L, {'Multiscale', 'Singlescale'}, 'Location', legendPos,...
                        'interpreter','latex', 'fontsize',15);
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'linear')
set(gca,'fontsize',15)

% subplot(2,1,2);
fig2 = figure;  
[dataFlip1] = AreaPlot(chagneSNR, lowerMSEsMulti', upperMSEsMulti');
[dataFlip2] = AreaPlot(chagneSNR, lowerMSEsSingle', upperMSEsSingle');

fill(dataFlip1.X, dataFlip1.Y, 'b', dataFlip2.X, dataFlip2.Y,'r','FaceAlpha',0.1);
hold on;
semilogy(chagneSNR, lowerMSEsMulti, 'b--',...
     chagneSNR, upperMSEsMulti, 'b--',...
     chagneSNR, upperMSEsSingle, 'r--',...
     chagneSNR, lowerMSEsSingle, 'r--')

semilogy(chagneSNR, meanMSEsMulti, 'b*-', chagneSNR, meanMSEsSingle, 'r*-');
ylabel('MSE', 'interpreter','latex')
% xlabel('p - Noise to Signal Ratio', 'interpreter','latex')
xlabel('SNR',  'interpreter','latex');
% xlim([factorss(1), factorss(end)]);
xlim([chagneSNR(end), chagneSNR(1)]);

if isSPDFlag
    ylim([minLim, maxLim]);
    yticks([0.005,0.01,0.02, 0.05])
    legendPos = 'Southwest';
    xticks([2, 5, 10, 40])
    xticklabels({'2', '5', '10', '40'})
else
    ylim([minLim, maxLim]);
    yticks([0.03,0.1,0.3])
    yticklabels({'0.03', '0.1', '0.3'})
    legendPos = 'Southwest';
    xticks([0.2,1,10])
    xticklabels({'0.2', '1', '10'})
end


L(1) = plot(nan, nan, 'b-');
L(2) = plot(nan, nan, 'r-');
legend(L, {'Multiscale', 'Singlescale'},'Location',legendPos, ...
           'interpreter','latex' ,'fontsize',15)

set(gca,'fontsize',15);
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')



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
