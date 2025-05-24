%%
clc;
close all;
%% Define test func & define imageres
if ~exist('imageRes', 'var')
   imageRes = 0;
end
      
fucnTest = AddGuassianNoiseAndOutliers(func, 0, 0, 0);

%% Initilize error strtucts
nLevels = length(ratesToSubSample);

errorsDelta1 = zeros(nRep, nLevels, 6);
errorsDelta2 = zeros(nRep, nLevels, 6);

nTestSample = (length(gridLimitsTest(1) : distTest : gridLimitsTest(2)))^nDim;

tradeoffDeltaMISE1 = zeros(nLevels, nTestSample, 3);
tradeoffDeltaMISE2 = zeros(nLevels, nTestSample, 3);

tradeoffInegrated1 = zeros(nLevels, 3);
tradeoffInegrated2 = zeros(nLevels, 3);

% meshNorms = zeros(nRep, nLevels); % Seems redundent

%% Create tmp data
bankOfEstimationsDetla1 = cell(nRep , 1); 
bankOfEstimationsDetla2 = cell(nRep , 1);
tmpMeshNorm = zeros(nRep, nLevels);

valueNormk = zeros(nRep, nLevels);

tmperalNoisedFunction = cell(nRep, 1);
dataSets = cell(nRep, 1);

%% Main loop
parfor iRep = 1 : nRep
    %% Add noise
    [funcNoised] = AddGuassianNoiseAndOutliers(func, sigma, pOut, sigmaOut);
%             [funcNoised] = AddOutliersAnat(func, pOut, hshOutlier, 2); % the last input is sceario number (1-3)
    tmperalNoisedFunction{iRep} = funcNoised;
    %% Sample source
    if flagScatter
        [pointsSet, samplingPlaces, outliersPos] = SampleNDfunctionScatterUniform(...
                                    funcNoised, gridLimits, dist, nDim);
    else
        [pointsSet, samplingPlaces, outliersPos] = SampleRnfunctionEqui(funcNoised,...
                gridLimits, dist, nDim);
    end
    dataSets{iRep} = pointsSet;
    %% Generate sub-Samples
    if flagScatter
        [samplingEachStage] = subSample2dScatterFirst(...
                    samplingPlaces, ratesToSubSample, nDim);
        forSingle = samplingEachStage;
%         [samplingEachStage] = subSample2dScatterEqualPartNoHeirarchical(...
%                 samplingPlaces, ratesToSubSample, nDim);
         
    else
        [samplingEachStage] = subSample2dRandom(samplingPlaces, ratesToSubSample);

    end
    %% Compute the mesh-nrom

    tmpMeshNorm(iRep,:) = ComputeMeshNormFromScatterDataRn(samplingPlaces,...
                            samplingEachStage, gridLimits, nDim);


    deltas = tmpMeshNorm(iRep,:) * nu;

    %% Multiscale
    [bankOfEstimationsDetla1{iRep},~, additionalInput] = MultiScale(...
             pointsSet, samplingEachStage, kernel, deltas, estFunc, nDim, spaceData);

                                              
    %% Single
    samplingEachStage1 = samplingEachStage;

    deltas1 = deltas;
    [bankOfEstimationsDetla2{iRep},~,~] = SingleScale(pointsSet,samplingEachStage1,...
                                                  kernel, deltas1, estFunc, nDim, spaceData);

    disp(['Finish: ' num2str(iRep)]);
end
disp('finish all Approximation')
% mean(tmpMeshNorm(:,3))
% median(tmpMeshNorm(:,3))

%% Initlaze to compute error


testFuncForAnalysis = repmat({fucnTest}, nRep, 1);

simplfiedAnalysisfunc = @(bankOfEstimations1, fucnTest1)...
ComputeErrorAccording3MeasuremntsByLevel(bankOfEstimations1, fucnTest1,...
                                  gridLimitsTest, distTest, nRep, nDim, spaceData);

%% compute for mutliscale without noise                              
[meanTmp, medianTmp, maxTmp] = simplfiedAnalysisfunc(...
                                bankOfEstimationsDetla1, testFuncForAnalysis);
errorsDelta1(:, :, 1) = meanTmp;
errorsDelta1(:, :, 2) = medianTmp;
errorsDelta1(:, :, 3) = maxTmp;

disp('Finish to compute error Multiscale')

%% compute for single without noise                              
[meanTmp, medianTmp, maxTmp] = simplfiedAnalysisfunc(...
                                bankOfEstimationsDetla2, testFuncForAnalysis);
errorsDelta2(:, :, 1) = meanTmp;
errorsDelta2(:, :, 2) = medianTmp;
errorsDelta2(:, :, 3) = maxTmp;

disp('Finish to compute error SingleScale')


%% Compute the acual error and not just measuremtns - Multiscale

if nDim == 1 % Currently 13.4.24, this information only used for grpahics 11.
% And it is only avilable in 1D. So this cut some running time.
    [errorMatMultiscale, pointTestForError] =...
        ComputeErrorAccording3MeasuremntsByPoint(bankOfEstimationsDetla1, fucnTest,...
                                          gridLimitsTest, distTest, nRep, nDim);

    [errorMatSinglescale, ~] =...
        ComputeErrorAccording3MeasuremntsByPoint(bankOfEstimationsDetla2, fucnTest,...
                                         gridLimitsTest, distTest, nRep, nDim);                                               
else
    errorMatMultiscale = [];
    errorMatSinglescale = [];
    pointTestForError = [];
end
                                      
%% Prepering for bias-variance tradeoff analysis
zeroFunc1 = @(x) 0 .* x;

%% update lambda
lambda = nu .* tmpMeshNorm;
%% Create Test Points 

nTestSample = (length(gridLimitsTest(1) : distTest : gridLimitsTest(2)))^nDim;

if nDim == 1
    testPoints = rand(nTestSample, nDim) .*...
                 (gridLimitsTest(2) - gridLimitsTest(1)) + gridLimitsTest(1);
    testPoints = sort(testPoints,1);
    error('not a good Factor');
elseif nDim == 2
    [XtestPoints, YtestPoints] = meshgrid(linspace(gridLimitsTest(1),...
                    gridLimitsTest(2), sqrt(nTestSample)));
    testPoints = [XtestPoints(:), YtestPoints(:)];

end
% Assume the inegrate dfcator is for a equidistance grid
intergrationFactor =  (gridLimitsTest(2) - gridLimitsTest(1)).^nDim ./ nTestSample;

%% MISE - Multiscale
% the output order is MSE, bias squared and variance (or the inegrated
% equivelnt)
[tradeoffDeltaMISE1, tradeoffInegrated1] = ComputeMSEMeusurmentes(...
            bankOfEstimationsDetla1, fucnTest, testPoints, nRep, intergrationFactor);

disp('Finish computing the Multiscale MSE')

%% MISE - Single Scale
[tradeoffDeltaMISE2, tradeoffInegrated2] = ComputeMSEMeusurmentes(...
            bankOfEstimationsDetla2, fucnTest, testPoints, nRep, intergrationFactor);

disp('Finish computing the Single MSE')
 
%% Create Snapshots

% Source
sourceSanpshoot = func(testPoints);

singleSanpshoot = zeros(size(sourceSanpshoot));
multiSanpshoot = zeros(size(sourceSanpshoot));

% Single & multi
funcSingleSnap = bankOfEstimationsDetla2{1};
funcSingleSnap = funcSingleSnap{end};

funcMultiSnap = bankOfEstimationsDetla1{1};
funcMultiSnap = funcMultiSnap{end};

parfor iPoint = 1 : size(testPoints,1)
    singleSanpshoot(iPoint) = funcSingleSnap(testPoints(iPoint,:));
    multiSanpshoot(iPoint) = funcMultiSnap(testPoints(iPoint,:));
end




%% disp and stuff

if sum(isnan(errorsDelta2(:))) || sum(isnan(errorsDelta1(:)))
    error('Problemmmm there is NAN')
end
clc;
%% plots error as function of mesh-nrom
% I assume that is used the same holton process in each
% meshs = squeeze(meshNorms(:, 1, 1, end));

% results = zeros(npOuts, nSigmas,3, 3);
% Ask user if she/he wants to save the figures
if ~toSave
    disp(['Hey, You are running my code. Your Saveing Path is ', savingPath]);
    disp('If you want to save the figure Press: 1 else 0');
    toSave = input('Enter 1 or 0: ');   
end

if(toSave)
        mkdir(savingPath);
        save([savingPath,'data.mat'], 'errorsDelta1', 'errorsDelta2',...
            'tradeoffDeltaMISE1', 'tradeoffDeltaMISE2',...
            'savingPath', 'nLevels',...
            'gridLimits','nRep', 'nDim', 'dist', 'mu', 'flagScatter',...
            'pointTestForError', 'errorMatMultiscale', 'errorMatSinglescale', 'testPoints',...
            'multiSanpshoot', 'sourceSanpshoot', 'singleSanpshoot', 'imageRes',...
            'tradeoffInegrated1', 'tradeoffInegrated2');
end

% X label
scale = [1];
for i = 1 : nLevels - 1
    scale = [mu .* scale(1) scale];
end
scale = scale * 100;
lebelsScale = 'Data [%]';


[figs, meanMulti, meanSingle] = Graphs1DEnvelopeNew(scale, lebelsScale,errorsDelta1(:,:,1:3), errorsDelta2(:,:,1:3));


if(toSave)
    if length(figs) == 1
        figs = {figs};
    end
    for i = 1 : length(figs)
        fileName = ['Evenlope_Ind_pOut_' num2str(i) 'RelError'];
        saveas(figs{i},[savingPath, fileName, '.fig']);
        saveas(figs{i},[savingPath, fileName, '.jpg']);

        fileName = ['Evenlope_Ind_pOut_' num2str(i) 'RelErrorVar'];
%         saveas(figsDiff{i},[savingPath, fileName, '.fig']);
%         saveas(figsDiff{i},[savingPath, fileName, '.jpg']);
    end
end
[figs] = Graphs10MSEByLevel1D(scale, lebelsScale, tradeoffDeltaMISE1, tradeoffDeltaMISE2,...
                         testPoints, nLevels, nDim, tradeoffInegrated1, tradeoffInegrated2...
                         , meanMulti, meanSingle);
if(toSave)
    if length(figs) == 1
        figs = {figs};
    end
    for i = 1 : length(figs)
        fileName = ['Evenlope_MSE_By_Level_' num2str(i)];
        saveas(figs{i},[savingPath, fileName, '.fig']);
        saveas(figs{i},[savingPath, fileName, '.jpg']);
        saveas(figs{i},[savingPath, fileName, '.pdf']);
    end
end


[figs] = Graphs11ErrorByPoint(pointTestForError, errorMatMultiscale,...
                            errorMatSinglescale, nLevels, nDim);
if(toSave)
    if length(figs) == 1
        figs = {figs};
    end
    for i = 1 : length(figs)
        fileName = ['Evenlope_Graph_Error_by_Error_Meusremets_' num2str(i)];
        saveas(figs{i},[savingPath, fileName, '.fig']);
        saveas(figs{i},[savingPath, fileName, '.jpg']);
    end
end


[figs] = Graphs12SnapShoot(nDim, testPoints,...
    multiSanpshoot, sourceSanpshoot, singleSanpshoot);
if(toSave)
    if length(figs) == 1
        figs = {figs};
    end
    for i = 1 : length(figs)
        fileName = ['SnapShoot_' num2str(i)];
        saveas(figs{i},[savingPath, fileName, '.fig']);
        saveas(figs{i},[savingPath, fileName, '.jpg']);
    end
end


