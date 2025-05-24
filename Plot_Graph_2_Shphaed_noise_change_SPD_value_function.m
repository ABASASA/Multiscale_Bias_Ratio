%%
clc;
close all;

addpath('kernels', 'EstSchemes', 'Sampling', ...
        'Auxiliary','MLS2D',  'GraphicalFunc',  'Error&MSE&MISE',...
        'NoiseingFunc', 'ComputeMetrics', 'SPD', 'Multiscale');

%% Envelope
sigmas = [0.05:0.05:0.25]; % In this run we doing by noise levels 
    %% Params

nDim = 2; % This parameter declar on which R^{nDim} we are working
nDimOut = 9; % The target function number dim

[spaceData] = SpaceDataSPD3();

gridLimits = [0, 1];  % define grid limies [a,b]^2
nRep = 100; % Repeat per source
dist = 0.1; % dist between sample (as it was in equidtiance grid, but it is scatter)
mu = 0.8; % dist between sample (as it was in equidtiance grid, but it is scatter)
nLevels = 3; % number of levels
ratesToSubSample = CreateRate(mu, nLevels, nDim); % Each time we take half of the field

distTest = 0.075; % the dist between grid of the test grid
gridLimitsTest = gridLimits + (distTest/2) * distTest + 2*[dist, -dist]; % This is not the distTest - this is anti-extrapolation

nu = 3 ; % For shephard

%% Oultieer params
pOut = 0; % probability
sigmaOut = 0; % Sigma
%% 
% In the following function path, you will save each combination results:
savingPathFunc = @(sigma) ['Graphs\Exp_2_sigma_' num2str(sigma) '_dist_' num2str(dist)...
                        '\'];

% In this path, you save the final plot:
savingPath = ['Graphs\Exp_2_dist_' num2str(dist) '_Unite\'];

                    
toSave = true;
flagScatter = true;

%% Kernel and estimation function
kernel = @(X1, X2, delta) Norm2KernelFunction(X1, X2, delta, nDim);
estFunc = @ShepperdSPDScehem;

%% Target Function
func = @(X) ExampleSPD(X);

%% Run
for i = 1 : length(sigmas)
    sigma = sigmas(i);
    savingPath = savingPathFunc(sigma);
    
    Envelope_Compare_Multiscale_Singlescale_SPD_value_function;
end    
%% Plot

factorss = sigmas;
funcPath = @(sigma) [savingPathFunc(sigma), 'data.mat'];

isSPDFlag = true; % This flag is just for minor grpahicks adjucjemts

Graph1ConbineNoises; 