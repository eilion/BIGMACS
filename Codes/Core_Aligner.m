function [results] = Core_Aligner(inputFile)

% Initialization:
path = ['Inputs/Core_Alignments/',inputFile,'/stack.txt'];
if exist(path,'file') ~= 2
    path = 'Defaults/stack.txt';
end
stack = getStack(path);

[data_type,IsLearn_transition,nSamples_learning,nSamples] = getSetting_Align(inputFile);

record_path = ['Inputs/Core_Alignments/',inputFile,'/records.txt'];

cal_curve = calibration_curve(stack(:,1));
data = getData(record_path,data_type);
[core_param,global_param] = Initializer(data,stack,cal_curve,['Inputs/Core_Alignments/',inputFile],data_type);
data = getAugData(core_param,data);
disp('## Parameters are initialized.');

[Samples,QQ] = init_Samples(data,core_param,global_param,cal_curve,stack,data_type);
disp('   Proposal distributions are initialized.');

% Iteration:
max_iters = 10;
threshold = 0.5;
iters = 0;
TT = ones(length(Samples),1);
while iters < max_iters && sum(TT) > 0
    iters = iters + 1;
    disp(['## Iteration ',num2str(iters),':']);
    
    [Samples,QQ] = get_Samples(Samples,QQ,data,core_param,global_param,cal_curve,stack,data_type,nSamples_learning,TT);
    disp('   Ages are sampled.');
    
    new_core_param = Update_Param_Transition(data,Samples,core_param,global_param,stack,IsLearn_transition);
    new_core_param = Update_Param_Shift_2(data,Samples,new_core_param,stack,data_type);
    
    loglik_old = getLOGLIK(data,Samples,global_param,core_param,cal_curve,stack,data_type);
    loglik_new = getLOGLIK(data,Samples,global_param,new_core_param,cal_curve,stack,data_type);
    
    new_core_param = Update_Param_R_2(data,Samples,global_param,new_core_param,stack,TT);
    
    disp('   Parameters are updated.');
    
    if strcmp(IsLearn_transition,'no') == 1 || strcmp(IsLearn_transition,'local') == 1
        diff = loglik_new - loglik_old;
        tt = ['   Updating parameters makes the average log-likelihood of samples increased by [',num2str(diff'),'].'];
    elseif strcmp(IsLearn_transition,'global') == 1
        diff = mean(loglik_new) - mean(loglik_old);
        tt = ['   Updating parameters makes the average log-likelihood of samples increased by ',num2str(diff),'.'];
    end
    disp(tt);
    
    core_param = new_core_param;
    
    if iters > 3
        if strcmp(IsLearn_transition,'no') == 1 || strcmp(IsLearn_transition,'local') == 1
            TT = (abs(loglik_old-loglik_new)>threshold);
        elseif strcmp(IsLearn_transition,'global') == 1
            if abs(diff) < threshold
                TT = zeros(length(Samples),1);
            end
        end
    end
end

% Sampling:
disp('## Parameters are learned. Sampling algorithm is now running...');
[Samples,~] = get_Samples(Samples,QQ,data,core_param,global_param,cal_curve,stack,data_type,nSamples,ones(length(Samples),1));
disp('   Ages are sampled.');
disp('----------------------------------------------------------------------------');


% Termination:
new_alignment_results_name = [inputFile,'_',data_type];
results = struct('name',cell(1,1),'type',cell(1,1),'data',cell(1,1),'samples',cell(1,1),'core_param',cell(1,1),'stack',cell(1,1));
results.name = new_alignment_results_name;
results.type = data_type;
results.data = data;
results.samples = Samples;
results.core_param = core_param;
results.stack = stack;
disp('## All procedures are done. The program is terminated.');

% Storing the results:
path = ['Outputs/Core_Alignments/',new_alignment_results_name];
if exist(path,'dir') == 7
    n = 0;
    DET = 1;
    while DET == 1
        n = n + 1;
        path = ['Outputs/Core_Alignments/',new_alignment_results_name,'(',num2str(n),')'];
        if exist(path,'dir') == 0
            DET = 0;
        end
    end
end
mkdir(path);

fileID = [path,'/results.mat'];
save(fileID,'results');

path_txt = [path,'/ages'];
mkdir(path_txt);
Info = getInfo(results);
for ll = 1:length(Info)
    fileID = [path_txt,'/',Info(ll).name,'.txt'];
    fid = fopen(fileID,'wt');
    fprintf(fid,'depth(m) lower_95(kyr) median(kyr) upper_95(kyr)');
    fprintf(fid,'\n');
    for i = 1:length(Info(ll).depth)
        fprintf(fid,'%f %f %f %f',[Info(ll).depth(i),Info(ll).lower(i),Info(ll).median(i),Info(ll).upper(i)]);
        fprintf(fid,'\n');
    end
    fclose(fid);
end

disp(['#  Results are stored in ',path,'.']);
disp('----------------------------------------------------------------------------');


end