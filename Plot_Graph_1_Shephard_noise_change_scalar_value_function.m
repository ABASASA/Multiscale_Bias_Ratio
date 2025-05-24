%%
clc;
close all;
clear;

addpath('kernels', 'EstSchemes', 'Sampling', ...
        'Auxiliary','MLS2D',  'GraphicalFunc',  'Error&MSE&MISE',...
        'NoiseingFunc', 'ComputeMetrics', 'Multiscale');
    
%% Envelope
sigmas = 0.25:0.25:2.5; % In this run we doing by noise levels 

    %% Params
nDim = 2 % This parameter declar on which R^{nDim} we are working
[spaceData] = SpaceDataEuclidian();
nDimOut = 1; % The target function number dim

gridLimits = [0, 1]; % define grid limies [a,b]^2
nRep = 100; % Repeat per source
dist = 0.075; % dist between sample (as it was in equidtiance grid, but it is scatter)
mu = 0.8; % The percent of same ddata betweeen scales of the multiscale
nLevels = 3; % number of levels
ratesToSubSample = CreateRate(mu, nLevels, nDim); % Each time we take half of the field

distTest = 0.02; % the dist between grid of the test grid
gridLimitsTest = gridLimits + (distTest/2) * distTest + 2*[dist, -dist]; % This is not the distTest - this is anti-extrapolation

nu = 3 ; % For shephard


%% Oultiers params - this part it optional and it nor used
pOut = 0; % probability - not used
sigmaOut = 0; % Sigma - not uesd

%% Saveing Path data

% In the following function path, you will save each combination results:
savingPathFunc = @(sigma) ['Graphs\Exp_1_sigma_' num2str(sigma) '_dist_' num2str(dist)...
                        '\'];

% In this path, you save the final plot:
savingPath = ['Graphs\Exp_1_dist_' num2str(dist) '_Unite\'];


toSave = true; % Flag - for automatic save
flagScatter = true; % notte if scattered data

%% Kernel and estimation function
kernel = @(X1, X2, delta) Norm2KernelFunction(X1, X2, delta, nDim);
estFunc = @ShepperdScehem;

%% Target Function

func =  @(X) sin(2 * pi * X(:,1)) ...
    .* (cos(4 * pi * X(:,2)));

%% Run
for i = 1 : length(sigmas)
    sigma = sigmas(i);
    savingPath = savingPathFunc(sigma);
    
    Envelope_Compare_Multiscale_Singlescale_scalar_value_function;
end    


%% Plot
factorss = sigmas;
funcPath = @(sigma) [savingPathFunc(sigma), 'data.mat'];
isSPDFlag = false; % This flag is just for minor grpahics adjucjemts

Graph1ConbineNoises;