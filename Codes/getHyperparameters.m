function [global_param] = getHyperparameters(inputFile_path)

global_param = struct('alpha',cell(1,1),'beta',cell(1,1),'q',cell(1,1),'d',cell(1,1),'nParticles',cell(1,1),'v',cell(1,1),'max_iters',cell(1,1),'h',cell(1,1),'tran_param',cell(1,1));

path = 'Defaults/hyperparameter.txt';
fileID = fopen(path);
INFO = textscan(fileID,'%s %s');
fclose(fileID);

global_param.alpha = str2double(INFO{2}(1));
global_param.beta = str2double(INFO{2}(2));
global_param.q = str2double(INFO{2}(3));
global_param.d = str2double(INFO{2}(4));
global_param.nParticles = str2double(INFO{2}(5));
global_param.v = str2double(INFO{2}(6));
global_param.max_iters = str2double(INFO{2}(7));
global_param.h = str2double(INFO{2}(8));

path = [inputFile_path,'/hyperparameter.txt'];
if exist(path,'file') == 2
    fileID = fopen(path);
    INFO = textscan(fileID,'%s %s');
    fclose(fileID);
    
    if sum(strcmp(INFO{1},'alpha:')==1) == 1
        global_param.alpha = str2double(INFO{2}{strcmp(INFO{1},'alpha:')==1});
    end
    
    if sum(strcmp(INFO{1},'beta:')==1) == 1
        global_param.beta = str2double(INFO{2}{strcmp(INFO{1},'beta:')==1});
    end
    
    if sum(strcmp(INFO{1},'q:')==1) == 1
        global_param.q = str2double(INFO{2}{strcmp(INFO{1},'q:')==1});
    end
    
    if sum(strcmp(INFO{1},'d:')==1) == 1
        global_param.d = str2double(INFO{2}{strcmp(INFO{1},'d:')==1});
    end
    
    if sum(strcmp(INFO{1},'nParticles:')==1) == 1
        global_param.nParticles = str2double(INFO{2}{strcmp(INFO{1},'nParticles:')==1});
    end
    
    if sum(strcmp(INFO{1},'v:')==1) == 1
        global_param.v = str2double(INFO{2}{strcmp(INFO{1},'v:')==1});
    end
    
    if sum(strcmp(INFO{1},'max_iters:')==1) == 1
        global_param.max_iters = str2double(INFO{2}{strcmp(INFO{1},'max_iters:')==1});
    end
    
    if sum(strcmp(INFO{1},'h:')==1) == 1
        global_param.h = str2double(INFO{2}{strcmp(INFO{1},'h:')==1});
    end
end


end