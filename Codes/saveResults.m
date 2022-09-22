function [savePath] = saveResults(data,samples,param,target,setting,inputFile,inputMode)

L = length(samples);

summary = struct('name',cell(L,1),'depth',cell(L,1));
for ll = 1:L
    summary(ll).name = data(ll).name;
    summary(ll).depth = data(ll).depth;
    summary(ll).d18O = data(ll).d18O;
    summary(ll).d18O_shift = data(ll).shift;
    summary(ll).d18O_scale = data(ll).scale;
    summary(ll).radiocarbon = data(ll).radiocarbon;
    summary(ll).additional_ages = data(ll).suggested_age;
    summary(ll).age_samples = samples(ll).ages;
    summary(ll).isoutlier = samples(ll).isoutlier;
    
    summary(ll).lower_95 = quantile(samples(ll).ages,0.025,2);
    summary(ll).lower_68 = quantile(samples(ll).ages,0.16,2);
    summary(ll).median = quantile(samples(ll).ages,0.5,2);
    summary(ll).mean = mean(samples(ll).ages,2);
    summary(ll).upper_68 = quantile(samples(ll).ages,0.84,2);
    summary(ll).upper_95 = quantile(samples(ll).ages,0.975,2);
    
    summary(ll).R = data(ll).R;
end


setting_alignment = struct('data_type',cell(1,1));
setting_alignment.data_type = setting.data_type;
setting_alignment.IsLearn_transition = setting.IsLearn_transition;
setting_alignment.stack_min = setting.stack_min;
setting_alignment.stack_max = setting.stack_max;


hyperparameter = struct('q',cell(1,1));
hyperparameter.q = param.q;
hyperparameter.d = param.d;
hyperparameter.nParticles = param.nParticles;
hyperparameter.max_iters = param.max_iters;
hyperparameter.a_d18O = param.a_d18O;
hyperparameter.b_d18O = param.b_d18O;
hyperparameter.a_C14 = param.a_C14;
hyperparameter.b_C14 = param.b_C14;
hyperparameter.nSamples_learning = setting.nSamples_learning;
hyperparameter.nSamples_drawing = setting.nSamples;


setting_core = struct('start_depth',cell(L,1));
for ll = 1:L
    path = 'Defaults/setting_core.txt';
    fileID = fopen(path);
    INFO = textscan(fileID,'%s %s');
    fclose(fileID);
    
    setting_core(ll).start_depth = str2double(INFO{2}{strcmp(INFO{1},'start_depth:')==1});
    setting_core(ll).end_depth = str2double(INFO{2}{strcmp(INFO{1},'end_depth:')==1});
    setting_core(ll).initial_shift = str2double(INFO{2}{strcmp(INFO{1},'initial_shift:')==1});
    setting_core(ll).initial_scale = str2double(INFO{2}{strcmp(INFO{1},'initial_scale:')==1});
    setting_core(ll).initial_average_sed_rate = str2double(INFO{2}{strcmp(INFO{1},'initial_average_sed_rate:')==1});
    setting_core(ll).islearn_shift = INFO{2}{strcmp(INFO{1},'islearn_shift:')==1};
    setting_core(ll).islearn_scale = INFO{2}{strcmp(INFO{1},'islearn_scale:')==1};
    setting_core(ll).islearn_average_sed_rate = INFO{2}{strcmp(INFO{1},'islearn_average_sed_rate:')==1};
    setting_core(ll).smoothness_bandwidth = str2double(INFO{2}{strcmp(INFO{1},'smoothness_bandwidth:')==1});
    setting_core(ll).lower_bound = str2double(INFO{2}{strcmp(INFO{1},'lower_bound:')==1});
    setting_core(ll).upper_bound = str2double(INFO{2}{strcmp(INFO{1},'upper_bound:')==1});
    setting_core(ll).min_resolution = str2double(INFO{2}{strcmp(INFO{1},'min_resolution:')==1});
    setting_core(ll).lower_sedrate = str2double(INFO{2}{strcmp(INFO{1},'lower_sedrate:')==1});
    setting_core(ll).upper_sedrate = str2double(INFO{2}{strcmp(INFO{1},'upper_sedrate:')==1});
    setting_core(ll).particle_bandwidth = str2double(INFO{2}{strcmp(INFO{1},'particle_bandwidth:')==1});
    
    path = ['Inputs/',inputFile,'/Records/',data(ll).name,'/setting_core.txt'];
    if exist(path,'file') == 2
        fileID = fopen(path);
        INFO = textscan(fileID,'%s %s');
        fclose(fileID);
        
        if sum(strcmp(INFO{1},'start_depth:')==1) == 1
            setting_core(ll).start_depth = str2double(INFO{2}{strcmp(INFO{1},'start_depth:')==1});
        end
        
        if sum(strcmp(INFO{1},'end_depth:')==1) == 1
            setting_core(ll).end_depth = str2double(INFO{2}{strcmp(INFO{1},'end_depth:')==1});
        end
        
        if sum(strcmp(INFO{1},'initial_shift:')==1) == 1
            setting_core(ll).initial_shift = str2double(INFO{2}{strcmp(INFO{1},'initial_shift:')==1});
        end
        
        if sum(strcmp(INFO{1},'initial_scale:')==1) == 1
            setting_core(ll).initial_scale = str2double(INFO{2}{strcmp(INFO{1},'initial_scale:')==1});
        end
        
        if sum(strcmp(INFO{1},'initial_average_sed_rate:')==1) == 1
            setting_core(ll).initial_average_sed_rate = str2double(INFO{2}{strcmp(INFO{1},'initial_average_sed_rate:')==1});
        end
        
        if sum(strcmp(INFO{1},'islearn_shift:')==1) == 1
            setting_core(ll).islearn_shift = INFO{2}{strcmp(INFO{1},'islearn_shift:')==1};
        end
        
        if sum(strcmp(INFO{1},'islearn_scale:')==1) == 1
            setting_core(ll).islearn_scale = INFO{2}{strcmp(INFO{1},'islearn_scale:')==1};
        end
        
        if sum(strcmp(INFO{1},'islearn_average_sed_rate:')==1) == 1
            setting_core(ll).islearn_average_sed_rate = INFO{2}{strcmp(INFO{1},'islearn_average_sed_rate:')==1};
        end
        
        if sum(strcmp(INFO{1},'lower_bound:')==1) == 1
            setting_core(ll).lower_bound = str2double(INFO{2}{strcmp(INFO{1},'lower_bound:')==1});
        end
        
        if sum(strcmp(INFO{1},'upper_bound:')==1) == 1
            setting_core(ll).upper_bound = str2double(INFO{2}{strcmp(INFO{1},'upper_bound:')==1});
        end
        
        if sum(strcmp(INFO{1},'min_resolution:')==1) == 1
            setting_core(ll).min_resolution = str2double(INFO{2}{strcmp(INFO{1},'min_resolution:')==1});
        end
        
        if sum(strcmp(INFO{1},'lower_sedrate:')==1) == 1
            setting_core(ll).lower_sedrate = str2double(INFO{2}{strcmp(INFO{1},'lower_sedrate:')==1});
        end
        
        if sum(strcmp(INFO{1},'upper_sedrate:')==1) == 1
            setting_core(ll).upper_sedrate = str2double(INFO{2}{strcmp(INFO{1},'upper_sedrate:')==1});
        end
        
        if sum(strcmp(INFO{1},'smoothness_bandwidth:')==1) == 1
            setting_core(ll).smoothness_bandwidth = str2double(INFO{2}{strcmp(INFO{1},'smoothness_bandwidth:')==1});
        end
        
        if sum(strcmp(INFO{1},'particle_bandwidth:')==1) == 1
            setting_core(ll).particle_bandwidth = str2double(INFO{2}{strcmp(INFO{1},'particle_bandwidth:')==1});
        end
    end
    
end



for ll = 1:L
    M = length(samples(ll).isoutlier);
    N = size(samples(ll).ages,1);
    D = size(data(ll).d18O,2);
    
    if strcmp(setting.data_type,'C14')
        AA = ones(N,M);
    else
        AA = zeros(N,M,D);
        for d = 1:D
            for m = 1:M
                AA(:,m,d) = samples(ll).isoutlier{m}(:,d);
            end
        end
    end
    samples(ll).isoutlier = AA;
    
    summary(ll).isoutlier = samples(ll).isoutlier;
end

if strcmp(inputMode,'alignment')
    
    data = rmfield(data,'is_top_14C_inlier');
    data = rmfield(data,'ACC_CONTRACTION');
    data = rmfield(data,'ACC_STEADY');
    data = rmfield(data,'ACC_EXPANSION');
    data = rmfield(data,'min_resolution_mode');
    data = rmfield(data,'res');
    data = rmfield(data,'min');
    data = rmfield(data,'max');
    data = rmfield(data,'SM_BW');
    data = rmfield(data,'PTCL_BW');
    
    TRAN_MAT = zeros(4,3);
    TRAN_MAT(1,:) = data(1).phi_I;
    TRAN_MAT(2,:) = data(1).phi_C;
    TRAN_MAT(3,:) = data(1).phi_M;
    TRAN_MAT(4,:) = data(1).phi_E;
    
    new_alignment_results_name = [inputFile,'_',setting.data_type];
    
    path = ['Outputs/',new_alignment_results_name];
    if exist(path,'dir') == 7
        n = 0;
        DET = 1;
        while DET == 1
            n = n + 1;
            path = ['Outputs/',new_alignment_results_name,'(',num2str(n),')'];
            if exist(path,'dir') == 0
                DET = 0;
            end
        end
    end
    mkdir(path);
    
    
    fileID = [path,'/transition_parameter.txt'];
    writematrix(TRAN_MAT,fileID,'Delimiter','tab');
    
    
    if strcmp(setting.data_type,'C14') || strcmp(setting.data_type,'both')
        CI_C14 = getCI_C14(data,target);
        
        path_txt = [path,'/C14_ages'];
        mkdir(path_txt);
        for ll = 1:length(CI_C14)
            fileID = [path_txt,'/',CI_C14(ll).name,'.txt'];
            fid = fopen(fileID,'wt');
            fprintf(fid,'depth(m) lower_95(kyr) lower_68(kyr) median(kyr) mean(kyr) upper_68(kyr) upper_95(kyr)');
            fprintf(fid,'\n');
            for i = 1:length(CI_C14(ll).depth)
                fprintf(fid,'%f %f %f %f %f %f %f',[CI_C14(ll).depth(i),CI_C14(ll).lower(i,:),CI_C14(ll).median(i,:),CI_C14(ll).upper(i,:)]);
                fprintf(fid,'\n');
            end
            fclose(fid);
        end
        
        
        for ll = 1:L
            CI_C14(ll).MD = CI_C14(ll).median;
        end
        CI_C14 = rmfield(CI_C14,'median');
        
        for ll = 1:L
            CI_C14(ll).lower_95 = CI_C14(ll).lower(:,1);
            CI_C14(ll).lower_68 = CI_C14(ll).lower(:,2);
            
            CI_C14(ll).median = CI_C14(ll).MD(:,1);
            CI_C14(ll).mean = CI_C14(ll).MD(:,2);
            
            CI_C14(ll).upper_68 = CI_C14(ll).upper(:,1);
            CI_C14(ll).upper_95 = CI_C14(ll).upper(:,2);
        end
        
        CI_C14 = rmfield(CI_C14,'lower');
        CI_C14 = rmfield(CI_C14,'MD');
        CI_C14 = rmfield(CI_C14,'upper');
        
        
        
        fileID = [path,'/results.mat'];
        save(fileID,'summary','CI_C14','target','setting_alignment','setting_core','hyperparameter');
    else
        fileID = [path,'/results.mat'];
        save(fileID,'summary','target','setting_alignment','setting_core','hyperparameter');
    end
    
    
    path_txt = [path,'/ages'];
    mkdir(path_txt);
    Info = getInfo(samples);
    for ll = 1:length(Info)
        fileID = [path_txt,'/',Info(ll).name,'.txt'];
        fid = fopen(fileID,'wt');
        fprintf(fid,'depth(m) lower_95(kyr) lower_68(kyr) median(kyr) mean(kyr) upper_68(kyr) upper_95(kyr)');
        fprintf(fid,'\n');
        for i = 1:length(Info(ll).depth)
            fprintf(fid,'%f %f %f %f %f %f %f',[Info(ll).depth(i),Info(ll).lower(i,:),Info(ll).median(i,:),Info(ll).upper(i,:)]);
            fprintf(fid,'\n');
        end
        fclose(fid);
    end
    
elseif strcmp(inputMode,'stacking')
    
    setting.start_age = setting.st;
    setting.end_age = setting.ed;
    
    setting_stacking = struct('start_age',cell(1,1));
    setting_stacking.start_age = setting.start_age;
    setting_stacking.end_age = setting.end_age;
    setting_stacking.interval = setting.interval;
    setting_stacking.interval_induced = setting.interval_induced;
    
    setting = rmfield(setting,'st');
    setting = rmfield(setting,'ed');
    
    param = rmfield(param,'eta');
    param = rmfield(param,'xi');
    param = rmfield(param,'lambda');
    param = rmfield(param,'lambda0');
    param = rmfield(param,'K');
    param = rmfield(param,'mu');
    param = rmfield(param,'nu');
    param = rmfield(param,'shift');
    param = rmfield(param,'scale');
    
    setting = rmfield(setting,'variance');
    setting = rmfield(setting,'kernel_function');
    
    data = rmfield(data,'is_top_14C_inlier');
    data = rmfield(data,'ACC_CONTRACTION');
    data = rmfield(data,'ACC_STEADY');
    data = rmfield(data,'ACC_EXPANSION');
    data = rmfield(data,'min_resolution_mode');
    data = rmfield(data,'res');
    data = rmfield(data,'min');
    data = rmfield(data,'max');
    data = rmfield(data,'SM_BW');
    data = rmfield(data,'PTCL_BW');
    
    TRAN_MAT = zeros(4,3);
    TRAN_MAT(1,:) = data(1).phi_I;
    TRAN_MAT(2,:) = data(1).phi_C;
    TRAN_MAT(3,:) = data(1).phi_M;
    TRAN_MAT(4,:) = data(1).phi_E;
    
    new_stack_name = [inputFile,'_',setting.data_type,'_stack'];
    
    path = ['Outputs/',new_stack_name];
    if exist(path,'dir') == 7
        n = 0;
        DET = 1;
        while DET == 1
            n = n + 1;
            path = ['Outputs/',new_stack_name,'(',num2str(n),')'];
            if exist(path,'dir') == 0
                DET = 0;
            end
        end
    end
    mkdir(path);
    
    fileID = [path,'/transition_parameter.txt'];
    writematrix(TRAN_MAT,fileID,'Delimiter','tab');
    
    if strcmp(setting.data_type,'C14') || strcmp(setting.data_type,'both')
        CI_C14 = getCI_C14(data,target);
        
        path_txt = [path,'/C14_ages'];
        mkdir(path_txt);
        for ll = 1:length(CI_C14)
            fileID = [path_txt,'/',CI_C14(ll).name,'.txt'];
            fid = fopen(fileID,'wt');
            fprintf(fid,'depth(m) lower_95(kyr) lower_68(kyr) median(kyr) mean(kyr) upper_68(kyr) upper_95(kyr)');
            fprintf(fid,'\n');
            for i = 1:length(CI_C14(ll).depth)
                fprintf(fid,'%f %f %f %f %f %f %f',[CI_C14(ll).depth(i),CI_C14(ll).lower(i,:),CI_C14(ll).median(i,:),CI_C14(ll).upper(i,:)]);
                fprintf(fid,'\n');
            end
            fclose(fid);
        end
        
        for ll = 1:L
            CI_C14(ll).MD = CI_C14(ll).median;
        end
        CI_C14 = rmfield(CI_C14,'median');
        
        for ll = 1:L
            CI_C14(ll).lower_95 = CI_C14(ll).lower(:,1);
            CI_C14(ll).lower_68 = CI_C14(ll).lower(:,2);
            
            CI_C14(ll).median = CI_C14(ll).MD(:,1);
            CI_C14(ll).mean = CI_C14(ll).MD(:,2);
            
            CI_C14(ll).upper_68 = CI_C14(ll).upper(:,1);
            CI_C14(ll).upper_95 = CI_C14(ll).upper(:,2);
        end
        
        CI_C14 = rmfield(CI_C14,'lower');
        CI_C14 = rmfield(CI_C14,'MD');
        CI_C14 = rmfield(CI_C14,'upper');
        
        
        MIN = inf;
        MAX = -inf;
        for ll = 1:L
            GG = quantile(summary(ll).age_samples(end,:),0.975);
            MAX = max(MAX,GG);
            
            GG = quantile(summary(ll).age_samples(1,:),0.025);
            MIN = min(MIN,GG);
        end
        
        MAX = ceil(MAX/10)*10;
        MIN = floor(MIN/10)*10;
        
        ID = (target.stack(:,1)<=MAX)&(target.stack(:,1)>=MIN);
        target.stack = target.stack(ID,:);
        target.stack_sample = target.stack_sample(ID,:,:);
        
        
        
        fileID = [path,'/results.mat'];
        save(fileID,'summary','CI_C14','target','setting_alignment','setting_stacking','setting_core','hyperparameter');
    else
        
        MIN = inf;
        MAX = -inf;
        for ll = 1:L
            GG = quantile(summary(ll).age_samples(end,:),0.975);
            MAX = max(MAX,GG);
            
            GG = quantile(summary(ll).age_samples(1,:),0.025);
            MIN = min(MIN,GG);
        end
        
        MAX = ceil(MAX/10)*10;
        MIN = floor(MIN/10)*10;
        
        ID = (target.stack(:,1)<=MAX)&(target.stack(:,1)>=MIN);
        target.stack = target.stack(ID,:);
        target.stack_sample = target.stack_sample(ID,:,:);
        
        
        fileID = [path,'/results.mat'];
        save(fileID,'summary','target','setting_alignment','setting_stacking','setting_core','hyperparameter');
    end
    
    
    fileID = [path,'/stack.txt'];
    fid = fopen(fileID,'wt');
    fprintf(fid,'age(kyr) mean(permil) sigma(permil)');
    fprintf(fid,'\n');
    for i = 1:size(target.stack,1)
        fprintf(fid,'%f %f %f',target.stack(i,:));
        fprintf(fid,'\n');
    end
    fclose(fid);
    
    fileID = [path,'/stack_samples_mean.txt'];
    writematrix(target.stack_sample(:,:,1),fileID,'Delimiter','space');
    
    fileID = [path,'/stack_samples_noisy.txt'];
    writematrix(target.stack_sample(:,:,2),fileID,'Delimiter','space');
    
    path_txt = [path,'/ages'];
    mkdir(path_txt);
    Info = getInfo(samples);
    for ll = 1:length(Info)
        fileID = [path_txt,'/',Info(ll).name,'.txt'];
        fid = fopen(fileID,'wt');
        fprintf(fid,'depth(m) lower_95(kyr) lower_68(kyr) median(kyr) mean(kyr) upper_68(kyr) upper_95(kyr)');
        fprintf(fid,'\n');
        for i = 1:length(Info(ll).depth)
            fprintf(fid,'%f %f %f %f %f %f %f',[Info(ll).depth(i),Info(ll).lower(i,:),Info(ll).median(i,:),Info(ll).upper(i,:)]);
            fprintf(fid,'\n');
        end
        fclose(fid);
    end
end

savePath = path;


end