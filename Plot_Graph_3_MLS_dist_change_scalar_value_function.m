%%
clc;
close all;
clear;

addpath('kernels', 'EstSchemes', 'Sampling', ...
        'Auxiliary','MLS2D',  'GraphicalFunc',  'Error&MSE&MISE',...
        'NoiseingFunc', 'ComputeMetrics', 'Multiscale');
%% Envelope
dists = round(logspace(log10(0.2),log10(0.025),7), 2);

%% Params
nDim = 2 % This parameter declar on which R^{nDim} we are working
[spaceData] = SpaceDataEuclidian();
nDimOut = 1; % The target function number dim

gridLimits = [0, 1]; % define grid limies [a,b]^2
nRep = 100; % Repeat per source

sigma = 0; % noise level
mu = 0.8; % dist between sample (as it was in equidtiance grid, but it is scatter)
nLevels = 3; % number of levels
ratesToSubSample = CreateRate(mu, nLevels, nDim); % Each time we take half of the field

distTest = 0.02; % the dist between grid of the test grid
gridLimitsTest = gridLimits + (distTest/2) * distTest + 2*[dist, -dist]; % This is not the distTest - this is anti-extrapolation

m = 1; % For moving least squares -- determent polynomial degree
nu = 3 ; 

%% Oultieer params
pOut = 0; % probability
sigmaOut = 0; % Sigma
%% 
% In the following function path, you will save each combination results:
savingPathFunc = @(dist) ['Graphs\Exp_3_dist_' num2str(dist) '\'];

% In this path, you save the final plot:
savingPath = ['Graphs\Exp_3_Unite\'];


toSave = true; % Flag - for automatic save
flagScatter = true;

%% Kernel and estimation function

kernel = @(X1, X2, delta) NormWeighForMLS2DExp(X1, X2, delta, nDim);%@Norm2WithRadLimitKernelFunction;
estFunc = @(XjMod, delta,...
                    kernel, additionalInput, nDim)MLS2DScheme(XjMod, delta,...
                    kernel, additionalInput, nDim, m);

%% Target Function

func = @(X) exp(sum(X.*(X), 2)) + 3; 

%% Run
for i = 1 : length(dists)
    dist = dists(i);
    savingPath = savingPathFunc(dist);
    
    Envelope_Compare_Multiscale_Singlescale_scalar_value_function;
end    
 

%% Plotting
funcPath =  @(dist) [savingPathFunc(dist), 'data.mat'];
factorss = dists;
xs =  (1 ./ dists).^2; % we use it to compute the number of points for each factor
isSPDFlag = false; % This flag is just for minor grpahicks adjucjemts

Graph1ConbineDist;
