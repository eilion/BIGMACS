function [data_type,kernel_function,IsLearn_transition,nSamples_learning,st,ed,interval] = getSetting_Stack(inputFile)

path = 'Defaults/setting_stack.txt';
fileID = fopen(path);
INFO = textscan(fileID,'%s %s');
fclose(fileID);

data_type = INFO{2}{1};
kernel_function = INFO{2}{2};
IsLearn_transition = INFO{2}{3};
nSamples_learning = str2double(INFO{2}{4});
st = str2double(INFO{2}{5});
ed = str2double(INFO{2}{6});
interval = str2double(INFO{2}{7});

path = ['Inputs/Stack_Construction/',inputFile,'/setting_stack.txt'];
if exist(path,'file') == 2
    fileID = fopen(path);
    INFO = textscan(fileID,'%s %s');
    fclose(fileID);
    
    M = 0;
    
    if sum(strcmp(INFO{1},'data_type:')==1) == 1
        data_type = INFO{2}{strcmp(INFO{1},'data_type:')==1};
        M = M + 1;
    end
    
    if sum(strcmp(INFO{1},'kernel_function:')==1) == 1
        kernel_function = INFO{2}{strcmp(INFO{1},'kernel_function:')==1};
        M = M + 1;
    end
    
    if sum(strcmp(INFO{1},'IsLearn_transition:')==1) == 1
        IsLearn_transition = INFO{2}{strcmp(INFO{1},'IsLearn_transition:')==1};
        M = M + 1;
    end
    
    if sum(strcmp(INFO{1},'nSamples:')==1) == 1
        nSamples_learning = str2double(INFO{2}{strcmp(INFO{1},'nSamples:')==1});
        M = M + 1;
    end
    
    if sum(strcmp(INFO{1},'start_age:')==1) == 1
        st = str2double(INFO{2}{strcmp(INFO{1},'start_age:')==1});
        M = M + 1;
    end
    
    if sum(strcmp(INFO{1},'end_age:')==1) == 1
        ed = str2double(INFO{2}{strcmp(INFO{1},'end_age:')==1});
        M = M + 1;
    end
    
    if sum(strcmp(INFO{1},'interval:')==1) == 1
        interval = str2double(INFO{2}{strcmp(INFO{1},'interval:')==1});
        M = M + 1;
    end
    
    if M < length(INFO{1})
        disp('   Warning: some stack settings in the input folder might not be specified as intended.');
    end
end


end