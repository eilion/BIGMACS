function BIGMACS(inputFile,inputMode)

disp('-------------------------------------------------------------------------------------------');
disp('## BIGMACS: Bayesian Inference Gaussian Process Multiproxy Alignment of Continuous Signals');
disp('-------------------------------------------------------------------------------------------');


if strcmp(inputMode,'alignment')
    disp('## Alignment algorithm is now running...');
    [data,samples,param,target,setting] = Alignment(inputFile);    
elseif strcmp(inputMode,'stacking')
    disp('## Stacking algorithm is now running...');
    [data,samples,param,target,setting] = Stacking(inputFile);
end


disp('#  Results and figures are being stored...');
savePath = saveResults(data,samples,param,target,setting,inputFile,inputMode);
saveFigures(inputFile,savePath);
disp('   Done.');
disp(['   Results and figures are stored in ',savePath,'.']);
disp('-------------------------------------------------------------------------------------------');


end