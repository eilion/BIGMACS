function [results] = Stack_Learner(inputFile)

% Initialization:
path = ['Inputs/Stack_Construction/',inputFile,'/initial_stack.txt'];
if exist(path,'file') ~= 2
    path = 'Defaults/stack.txt';
end
stack = getStack(path);

[data_type,kernel_function,IsLearn_transition,nSamples_learning,st,ed,interval] = getSetting_Stack(inputFile);

record_path = ['Inputs/Stack_Construction/',inputFile,'/records.txt'];

st = max([st,stack(1,1)]);
ed = min([ed,stack(end,1)]);
age = (st:interval:ed)';
stack = init_Stack(stack,age);
initial_stack = stack;
cal_curve = calibration_curve(stack(:,1));
data = getData(record_path,data_type);
[core_param,global_param] = Initializer(data,stack,cal_curve,['Inputs/Stack_Construction/',inputFile],data_type);
data = getAugData(core_param,data);
disp('## Parameters are initialized.');

[Samples,QQ] = init_Samples(data,core_param,global_param,cal_curve,stack,data_type);
disp('   Proposal distributions are initialized.');


% Iteration:
max_iters = 5;
max_iters2 = 10;
threshold = 0.5;
iters = 0;
while iters < max_iters
    iters = iters + 1;
    tt = ['## Iteration ',num2str(iters),':'];
    disp(tt);
    
    iters2 = 0;
    TT = ones(length(Samples),1);
    while iters2 < max_iters2 && sum(TT) > 0
        iters2 = iters2 + 1;
        
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
        
        if iters == 1
            if iters2 > 3
                if strcmp(IsLearn_transition,'no') == 1 || strcmp(IsLearn_transition,'local') == 1
                    TT = (abs(loglik_old-loglik_new)>threshold);
                elseif strcmp(IsLearn_transition,'global') == 1
                    if abs(diff) < threshold
                        TT = zeros(length(Samples),1);
                    end
                end
            end
        else
            if strcmp(IsLearn_transition,'no') == 1 || strcmp(IsLearn_transition,'local') == 1
                TT = (abs(loglik_old-loglik_new)>threshold);
            elseif strcmp(IsLearn_transition,'global') == 1
                if abs(diff) < threshold
                    TT = zeros(length(Samples),1);
                end
            end
        end
    end
    disp('   Parameters are learned.');
    [stack,rcd_stack] = Stack_Updater_2(Samples,data,initial_stack,stack,core_param,kernel_function);
    disp('   The stack is updated.');
end
disp('----------------------------------------------------------------------------');

% Termination:
new_stack_name = [inputFile,'_',data_type,'_',kernel_function];
results = struct('name',cell(1,1),'type',cell(1,1),'data',cell(1,1),'samples',cell(1,1),'core_param',cell(1,1),'rcd_stack',cell(1,1),'old_stack',cell(1,1),'stack',cell(1,1));
results.name = new_stack_name;
results.type = data_type;
results.data = data;
results.samples = Samples;
results.core_param = core_param;
results.rcd_stack = rcd_stack;
results.old_stack = initial_stack;
results.stack = stack;
disp('## All procedures are done. The program is terminated.');


% Storing the results:
path = ['Outputs/Stack_Construction/',new_stack_name];
if exist(path,'dir') == 7
    n = 0;
    DET = 1;
    while DET == 1
        n = n + 1;
        path = ['Outputs/Stack_Construction/',new_stack_name,'(',num2str(n),')'];
        if exist(path,'dir') == 0
            DET = 0;
        end
    end
end
mkdir(path);

fileID = [path,'/results.mat'];
save(fileID,'results');

fileID = [path,'/stack.txt'];
fid = fopen(fileID,'wt');
fprintf(fid,'age(kyr) mean(permil) sigma(permil)');
fprintf(fid,'\n');
for i = 1:size(stack,1)
    fprintf(fid,'%f %f %f',stack(i,:));
    fprintf(fid,'\n');
end
fclose(fid);

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