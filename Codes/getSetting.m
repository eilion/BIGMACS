function [setting] = getSetting(inputFile,MODE)

setting = struct('data_type',cell(1,1));

path = 'Defaults/setting_alignment.txt';
fileID = fopen(path);
INFO = textscan(fileID,'%s %s');
fclose(fileID);

setting.data_type = INFO{2}{strcmp(INFO{1},'data_type:')==1};
setting.IsLearn_transition = INFO{2}{strcmp(INFO{1},'IsLearn_transition:')==1};
setting.stack_min = str2double(INFO{2}{strcmp(INFO{1},'stack_min:')==1});
setting.stack_max = str2double(INFO{2}{strcmp(INFO{1},'stack_max:')==1});

path = 'Defaults/hyperparameter.txt';
fileID = fopen(path);
INFO = textscan(fileID,'%s %s');
fclose(fileID);
setting.nSamples_learning = str2double(INFO{2}{strcmp(INFO{1},'nSamples_learning:')==1});
setting.nSamples = str2double(INFO{2}{strcmp(INFO{1},'nSamples_drawing:')==1});


path = ['Inputs/',inputFile,'/setting_alignment.txt'];
if exist(path,'file') == 2
    fileID = fopen(path);
    INFO = textscan(fileID,'%s %s');
    fclose(fileID);
    
    if sum(strcmp(INFO{1},'data_type:')==1) == 1
        setting.data_type = INFO{2}{strcmp(INFO{1},'data_type:')==1};
    end
    
    if sum(strcmp(INFO{1},'IsLearn_transition:')==1) == 1
        setting.IsLearn_transition = INFO{2}{strcmp(INFO{1},'IsLearn_transition:')==1};
    end
    
    if sum(strcmp(INFO{1},'stack_min:')==1) == 1
        setting.stack_min = str2double(INFO{2}{strcmp(INFO{1},'stack_min:')==1});
    end
    
    if sum(strcmp(INFO{1},'stack_max:')==1) == 1
        setting.stack_max = str2double(INFO{2}{strcmp(INFO{1},'stack_max:')==1});
    end
end

path = ['Inputs/',inputFile,'/hyperparameter.txt'];
if exist(path,'file') == 2
    fileID = fopen(path);
    INFO = textscan(fileID,'%s %s');
    fclose(fileID);
    
    if sum(strcmp(INFO{1},'nSamples_learning:')==1) == 1
        setting.nSamples_learning = str2double(INFO{2}{strcmp(INFO{1},'nSamples_learning:')==1});
    end
    
    if sum(strcmp(INFO{1},'nSamples_drawing:')==1) == 1
        setting.nSamples = str2double(INFO{2}{strcmp(INFO{1},'nSamples_drawing:')==1});
    end
end


if strcmp(MODE,'stacking')
    path = 'Defaults/setting_stacking.txt';
    fileID = fopen(path);
    INFO = textscan(fileID,'%s %s');
    fclose(fileID);
    
    setting.variance = 'heteroscedastic';
    setting.kernel_function = 'OU';
    setting.st = str2double(INFO{2}{strcmp(INFO{1},'start_age:')==1});
    setting.ed = str2double(INFO{2}{strcmp(INFO{1},'end_age:')==1});
    setting.interval = str2double(INFO{2}{strcmp(INFO{1},'interval:')==1});
    setting.interval_induced = 0.5;
    
    path = ['Inputs/',inputFile,'/setting_stacking.txt'];
    if exist(path,'file') == 2
        fileID = fopen(path);
        INFO = textscan(fileID,'%s %s');
        fclose(fileID);
        
        if sum(strcmp(INFO{1},'start_age:')==1) == 1
            setting.st = str2double(INFO{2}{strcmp(INFO{1},'start_age:')==1});
        end
        
        if sum(strcmp(INFO{1},'end_age:')==1) == 1
            setting.ed = str2double(INFO{2}{strcmp(INFO{1},'end_age:')==1});
        end
        
        if sum(strcmp(INFO{1},'interval:')==1) == 1
            setting.interval = str2double(INFO{2}{strcmp(INFO{1},'interval:')==1});
        end

        if sum(strcmp(INFO{1},'interval_induced:')==1) == 1
            setting.interval_induced = str2double(INFO{2}{strcmp(INFO{1},'interval_induced:')==1});
        end
    end
end


end