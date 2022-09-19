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
end

if strcmp(inputMode,'alignment')
    
    setting = rmfield(setting,'nSamples_learning');
    setting = rmfield(setting,'nSamples');
    
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
    
    
    fileID = [path,'/new_transition_parameter.txt'];
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
        save(fileID,'summary','data','samples','CI_C14','param','target','setting');
    else
        fileID = [path,'/results.mat'];
        save(fileID,'summary','data','samples','param','target','setting');
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
    
    setting = rmfield(setting,'nSamples_learning');
    setting = rmfield(setting,'nSamples');
    
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
    
    fileID = [path,'/new_transition_parameter.txt'];
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
        save(fileID,'summary','data','samples','CI_C14','param','target','setting');
    else
        fileID = [path,'/results.mat'];
        save(fileID,'summary','data','samples','param','target','setting');
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