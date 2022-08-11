function [data,samples,param,target,setting] = Alignment(inputFile)

disp('#  Initializing...');
setting = getSetting(inputFile,'alignment');
target = getInitialTarget(inputFile,setting,'alignment');
[data,data_part,param,setting] = getData(inputFile,target,setting,'alignment');


% initialize samples:
QQ = initializeAlignment(data,data_part,param,target,setting);
disp('   Done.');


% alignment:
disp('#  Records are being aligned to the target stack...');
[data,samples,param] = getAlignment(data,data_part,QQ,param,target,setting,'alignment');
disp('   Done.');
disp('-------------------------------------------------------------------------------------------');


end