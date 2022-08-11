function [savePath] = saveResults(data,samples,param,target,setting,inputFile,inputMode)

if strcmp(inputMode,'alignment')
    
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
    
    
    if strcmp(setting.data_type,'C14') || strcmp(setting.data_type,'both')
        CI_C14 = getCI_C14(data);
        
        path_txt = [path,'/C14_ages'];
        mkdir(path_txt);
        for ll = 1:length(CI_C14)
            fileID = [path_txt,'/',CI_C14(ll).name,'.txt'];
            fid = fopen(fileID,'wt');
            fprintf(fid,'depth(m) lower_95.4(kyr) lower_68.2(kyr) median(kyr) mean(kyr) upper_68.2(kyr) upper_95.4(kyr)');
            fprintf(fid,'\n');
            for i = 1:length(CI_C14(ll).depth)
                fprintf(fid,'%f %f %f %f %f %f %f',[CI_C14(ll).depth(i),CI_C14(ll).lower(i,:),CI_C14(ll).median(i,:),CI_C14(ll).upper(i,:)]);
                fprintf(fid,'\n');
            end
            fclose(fid);
        end
        
        fileID = [path,'/results.mat'];
        save(fileID,'data','samples','CI_C14','param','target','setting');
    else
        fileID = [path,'/results.mat'];
        save(fileID,'data','samples','param','target','setting');
    end
    
    
    path_txt = [path,'/ages'];
    mkdir(path_txt);
    Info = getInfo(samples);
    for ll = 1:length(Info)
        fileID = [path_txt,'/',Info(ll).name,'.txt'];
        fid = fopen(fileID,'wt');
        fprintf(fid,'depth(m) lower_95.4(kyr) lower_68.2(kyr) median(kyr) mean(kyr) upper_68.2(kyr) upper_95.4(kyr)');
        fprintf(fid,'\n');
        for i = 1:length(Info(ll).depth)
            fprintf(fid,'%f %f %f %f %f %f %f',[Info(ll).depth(i),Info(ll).lower(i,:),Info(ll).median(i,:),Info(ll).upper(i,:)]);
            fprintf(fid,'\n');
        end
        fclose(fid);
    end

elseif strcmp(inputMode,'stacking')
    
    new_stack_name = [inputFile,'_',setting.data_type,'_',setting.kernel_function,'_',setting.variance];
    
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
    
    if strcmp(setting.data_type,'C14') || strcmp(setting.data_type,'both')
        CI_C14 = getCI_C14(data);
        
        path_txt = [path,'/C14_ages'];
        mkdir(path_txt);
        for ll = 1:length(CI_C14)
            fileID = [path_txt,'/',CI_C14(ll).name,'.txt'];
            fid = fopen(fileID,'wt');
            fprintf(fid,'depth(m) lower_95.4(kyr) lower_68.2(kyr) median(kyr) mean(kyr) upper_68.2(kyr) upper_95.4(kyr)');
            fprintf(fid,'\n');
            for i = 1:length(CI_C14(ll).depth)
                fprintf(fid,'%f %f %f %f %f %f %f',[CI_C14(ll).depth(i),CI_C14(ll).lower(i,:),CI_C14(ll).median(i,:),CI_C14(ll).upper(i,:)]);
                fprintf(fid,'\n');
            end
            fclose(fid);
        end
        
        fileID = [path,'/results.mat'];
        save(fileID,'data','samples','CI_C14','param','target','setting');
    else
        fileID = [path,'/results.mat'];
        save(fileID,'data','samples','param','target','setting');
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
        fprintf(fid,'depth(m) lower_95.4(kyr) lower_68.2(kyr) median(kyr) mean(kyr) upper_68.2(kyr) upper_95.4(kyr)');
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