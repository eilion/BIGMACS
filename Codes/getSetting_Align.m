function [data_type,IsLearn_transition,nSamples_learning,nSamples] = getSetting_Align(inputFile)

path = 'Defaults/setting_align.txt';
fileID = fopen(path);
INFO = textscan(fileID,'%s %s');
fclose(fileID);

data_type = INFO{2}{1};
IsLearn_transition = INFO{2}{2};
nSamples_learning = str2double(INFO{2}{3});
nSamples = str2double(INFO{2}{4});

path = ['Inputs/Core_Alignments/',inputFile,'/setting_align.txt'];
if exist(path,'file') == 2
    fileID = fopen(path);
    INFO = textscan(fileID,'%s %s');
    fclose(fileID);
    
    M = 0;
    
    if sum(strcmp(INFO{1},'data_type:')==1) == 1
        data_type = INFO{2}{strcmp(INFO{1},'data_type:')==1};
        M = M + 1;
    end
    
    if sum(strcmp(INFO{1},'IsLearn_transition:')==1) == 1
        IsLearn_transition = INFO{2}{strcmp(INFO{1},'IsLearn_transition:')==1};
        M = M + 1;
    end
    
    if sum(strcmp(INFO{1},'nSamples_learning:')==1) == 1
        nSamples_learning = str2double(INFO{2}{strcmp(INFO{1},'nSamples_learning:')==1});
        M = M + 1;
    end
    
    if sum(strcmp(INFO{1},'nSamples_drawing:')==1) == 1
        nSamples = str2double(INFO{2}{strcmp(INFO{1},'nSamples_drawing:')==1});
        M = M + 1;
    end
    
    if M < length(INFO{1})
        disp('   Warning: some alignment settings in the input folder might not be specified as intended.');
    end
end


end