function BIGMACS(inputFile,inputMode)

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
saveFigures(inputFile,savePath);
disp('   Done.');
disp(['   Results and figures are stored in ',savePath,'.']);
disp('-------------------------------------------------------------------------------------------');


end