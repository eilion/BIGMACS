function [initial_shift,initial_scale,initial_average_sed_rate,islearn_shift,islearn_scale,islearn_average_sed_rate,lower_bound,upper_bound,min_resolution,lower_sedrate,upper_sedrate] = getSetting_Core(core_name)

path = 'Defaults/setting_core.txt';
fileID = fopen(path);
INFO = textscan(fileID,'%s %s');
fclose(fileID);

initial_shift = str2double(INFO{2}{1});
initial_scale = str2double(INFO{2}{2});
initial_average_sed_rate = str2double(INFO{2}{3});
islearn_shift = INFO{2}{4};
islearn_scale = INFO{2}{5};
islearn_average_sed_rate = INFO{2}{6};
lower_bound = str2double(INFO{2}{7});
upper_bound = str2double(INFO{2}{8});
min_resolution = str2double(INFO{2}{9});
lower_sedrate = str2double(INFO{2}{10});
upper_sedrate = str2double(INFO{2}{11});


path = ['Cores/',core_name,'/setting_core.txt'];
if exist(path,'file') == 2
    fileID = fopen(path);
    INFO = textscan(fileID,'%s %s');
    fclose(fileID);
    
    M = 0;
    
    if sum(strcmp(INFO{1},'initial_shift:')==1) == 1
        initial_shift = str2double(INFO{2}{strcmp(INFO{1},'initial_shift:')==1});
        M = M + 1;
    end
    
    if sum(strcmp(INFO{1},'initial_scale:')==1) == 1
        initial_scale = str2double(INFO{2}{strcmp(INFO{1},'initial_scale:')==1});
        M = M + 1;
    end
    
    if sum(strcmp(INFO{1},'initial_average_sed_rate:')==1) == 1
        initial_average_sed_rate = str2double(INFO{2}{strcmp(INFO{1},'initial_average_sed_rate:')==1});
        M = M + 1;
    end
    
    if sum(strcmp(INFO{1},'islearn_shift:')==1) == 1
        islearn_shift = INFO{2}{strcmp(INFO{1},'islearn_shift:')==1};
        M = M + 1;
    end
    
    if sum(strcmp(INFO{1},'islearn_scale:')==1) == 1
        islearn_scale = INFO{2}{strcmp(INFO{1},'islearn_scale:')==1};
        M = M + 1;
    end
    
    if sum(strcmp(INFO{1},'islearn_average_sed_rate:')==1) == 1
        islearn_average_sed_rate = INFO{2}{strcmp(INFO{1},'islearn_average_sed_rate:')==1};
        M = M + 1;
    end
    
    if sum(strcmp(INFO{1},'lower_bound:')==1) == 1
        lower_bound = str2double(INFO{2}{strcmp(INFO{1},'lower_bound:')==1});
        M = M + 1;
    end
    
    if sum(strcmp(INFO{1},'upper_bound:')==1) == 1
        upper_bound = str2double(INFO{2}{strcmp(INFO{1},'upper_bound:')==1});
        M = M + 1;
    end
    
    if sum(strcmp(INFO{1},'min_resolution:')==1) == 1
        min_resolution = str2double(INFO{2}{strcmp(INFO{1},'min_resolution:')==1});
        M = M + 1;
    end
    
    if sum(strcmp(INFO{1},'lower_sedrate:')==1) == 1
        lower_sedrate = str2double(INFO{2}{strcmp(INFO{1},'lower_sedrate:')==1});
        M = M + 1;
    end
    
    if sum(strcmp(INFO{1},'upper_sedrate:')==1) == 1
        upper_sedrate = str2double(INFO{2}{strcmp(INFO{1},'upper_sedrate:')==1});
        M = M + 1;
    end
    
    if M < length(INFO{1})
        disp('   Warning: some core settings in the input folder might not be specified as intended.');
    end
end


end