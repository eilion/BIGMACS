function BIGMACS(varargin)
% function BIGMACS(inputFile,inputMode)

inputFile = varargin{1};
inputMode = varargin{2};

if length(varargin) == 2
    FIG_MODE = 'no';
elseif length(varargin) == 3
    FIG_MODE = varargin{3};
end


addpath('Codes/');

disp('-------------------------------------------------------------------------------------------');
disp('## BIGMACS: Bayesian Inference Gaussian Process Multiproxy Alignment of Continuous Signals');
disp('-------------------------------------------------------------------------------------------');


if strcmp(inputMode,'age_model_construction')
    inputMode = 'alignment';
    disp('## Age model construction algorithm is now running...');
    [data,samples,param,target,setting] = Alignment(inputFile);    
elseif strcmp(inputMode,'stack_construction')
    inputMode = 'stacking';
    disp('## Stack construction algorithm is now running...');
    [data,samples,param,target,setting] = Stacking(inputFile);
end


disp('#  Results and figures are being stored...');
savePath = saveResults(data,samples,param,target,setting,inputFile,inputMode);
saveFigures(inputFile,savePath,inputMode,FIG_MODE);
disp('   Done.');
disp(['   Results and figures are stored in ',savePath,'.']);
disp('-------------------------------------------------------------------------------------------');


end