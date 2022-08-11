function [data,samples,param,target,setting] = Stacking(inputFile)

disp('#  Initializing...');
setting = getSetting(inputFile,'stacking');
target = getInitialTarget(inputFile,setting,'stacking');
[data,data_part,param,setting] = getData(inputFile,target,setting,'stacking');


% initialize samples:
QQ = initializeAlignment(data,data_part,param,target,setting);
disp('   Done.');


% stacking:
disp('#  Stacking algorithm is now running...');
for iters = 1:5
    disp(['#  Iteration ',num2str(iters),':']);
    [data,samples,param] = getAlignment(data,data_part,QQ,param,target,setting,'stacking');
    disp('   Records are aligned to the target stack.');
    [param,target] = updateStack(data,samples,param,target,setting);
    disp('   A new stack is updated.');
    [param,target] = standardizeStack(data,param,target);
end
disp('-------------------------------------------------------------------------------------------');


end