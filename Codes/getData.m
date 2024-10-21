function [data,data_ps,param,setting] = getData(inputFile,target,setting,MODE)

% global parameters:
param = struct('q',cell(1,1),'d',cell(1,1),'nParticles',cell(1,1),'max_iters',cell(1,1),'tran_param',cell(1,1));
path = 'Defaults/hyperparameter.txt';
fileID = fopen(path);
INFO = textscan(fileID,'%s %s');
fclose(fileID);

param.q = str2double(INFO{2}{strcmp(INFO{1},'q:')==1});
param.d = str2double(INFO{2}{strcmp(INFO{1},'d:')==1});
param.nParticles = str2double(INFO{2}{strcmp(INFO{1},'nParticles:')==1});
param.max_iters = str2double(INFO{2}{strcmp(INFO{1},'max_iters:')==1});
param.a_d18O = str2double(INFO{2}{strcmp(INFO{1},'a_d18O:')==1});
param.b_d18O = str2double(INFO{2}{strcmp(INFO{1},'b_d18O:')==1});
param.a_C14 = str2double(INFO{2}{strcmp(INFO{1},'a_C14:')==1});
param.b_C14 = str2double(INFO{2}{strcmp(INFO{1},'b_C14:')==1});

path = ['Inputs/',inputFile,'/hyperparameter.txt'];
if exist(path,'file') == 2
    fileID = fopen(path);
    INFO = textscan(fileID,'%s %s');
    fclose(fileID);
    
    if sum(strcmp(INFO{1},'q:')==1) == 1
        param.q = str2double(INFO{2}{strcmp(INFO{1},'q:')==1});
    end
    
    if sum(strcmp(INFO{1},'d:')==1) == 1
        param.d = str2double(INFO{2}{strcmp(INFO{1},'d:')==1});
    end
    
    if sum(strcmp(INFO{1},'nParticles:')==1) == 1
        param.nParticles = str2double(INFO{2}{strcmp(INFO{1},'nParticles:')==1});
    end
    
    if sum(strcmp(INFO{1},'max_iters:')==1) == 1
        param.max_iters = str2double(INFO{2}{strcmp(INFO{1},'max_iters:')==1});
    end
    
    if sum(strcmp(INFO{1},'a_d18O:')==1) == 1
        param.a_d18O = str2double(INFO{2}{strcmp(INFO{1},'a_d18O:')==1});
    end
    
    if sum(strcmp(INFO{1},'b_d18O:')==1) == 1
        param.b_d18O = str2double(INFO{2}{strcmp(INFO{1},'b_d18O:')==1});
    end
    
    if sum(strcmp(INFO{1},'a_C14:')==1) == 1
        param.a_C14 = str2double(INFO{2}{strcmp(INFO{1},'a_C14:')==1});
    end
    
    if sum(strcmp(INFO{1},'b_C14:')==1) == 1
        param.b_C14 = str2double(INFO{2}{strcmp(INFO{1},'b_C14:')==1});
    end
end

if strcmp(MODE,'stacking')
    param.eta = 1;
    param.xi = sqrt(0.5/setting.interval_induced);
    param.lambda = 1;
    param.lambda0 = 1e-3;
    param.K = 10*ones(setting.nSamples_learning,1);
    param.mu = cell(setting.nSamples_learning,1);
    param.nu = cell(setting.nSamples_learning,1);
end

path = ['Inputs/',inputFile,'/transition_parameter.txt'];
if exist(path,'file') ~= 2
    path = 'Defaults/transition_parameter.txt';
end
tran_param = load(path);

param.tran_param = tran_param;
param.shift = 0;
param.scale = 1;


% Records:
path = ['Inputs/',inputFile,'/records.txt'];
if exist(path,'file') == 2
    core_title = textread(path,'%s');
    L = length(core_title);
    
    path = ['Inputs/',inputFile,'/Records'];
    list = dir(path);
    
    ID = zeros(L,1);
    for ll = 1:L
        if sum(strcmp({list.name},core_title{ll})) == 0
            disp(['   Warning: record ',core_title{ll},' is not found in Inputs/',inputFile,'/Records/.']);
            ID(ll) = 1;
        end
    end
    core_title(ID==1) = [];
    L = length(core_title);
else
    path = ['Inputs/',inputFile,'/Records'];
    list = dir(path);
    list(strcmp({list.name},'.')) = [];
    list(strcmp({list.name},'..')) = [];
    list(strcmp({list.name},'.DS_Store')) = [];
    
    L = length(list);
    core_title = cell(L,1);
    for ll = 1:L
        core_title{ll} = list(ll).name;
    end
end

data = struct('name',cell(L,1),'depth',cell(L,1),'d18O',cell(L,1),'radiocarbon',cell(L,1),'suggested_age',cell(L,1),'initial_age',cell(L,1),'scale',cell(L,1),'shift',cell(L,1),'R',cell(L,1),'phi_I',cell(L,1),'phi_C',cell(L,1),'phi_M',cell(L,1),'phi_E',cell(L,1),'min',cell(L,1),'max',cell(L,1),'res',cell(L,1),'lower_sedrate',cell(L,1),'upper_sedrate',cell(L,1));

for ll = 1:L
    data(ll).name = core_title{ll};
    
    % Gather d18O data:
    path = ['Inputs/',inputFile,'/Records/',data(ll).name,'/d18O_data.txt'];
    if strcmp(setting.data_type,'C14') == 0 && exist(path,'file') == 2
        fileID = fopen(path);
        AA = textscan(fileID,'%s %s');
        fclose(fileID);
        d18O_data_temp = zeros(length(AA{1})-1,2);
        for k = 1:2
            d18O_data_temp(:,k) = str2double(AA{k}(2:end));
        end
        [~,order] = sort(d18O_data_temp(:,1),'ascend');
        d18O_data_temp = d18O_data_temp(order,:);
        depth_d18O_temp = d18O_data_temp(:,1);
        d18O_data_temp = d18O_data_temp(:,2);
        
        N = length(depth_d18O_temp);
        MAX = zeros(N,1);
        for n = 1:N
            MAX(n) = sum(abs(depth_d18O_temp(n)-depth_d18O_temp)<1e-8);
        end
        MAX = max(MAX);
        
        d18O_depth = zeros(N,1);
        d18O_data = NaN*ones(N,MAX);
        
        count = 1;
        d18O_depth(count) = depth_d18O_temp(1);
        d18O_data(count,1) = d18O_data_temp(1);
        count2 = 1;
        for n = 2:N
            if abs(depth_d18O_temp(n)-d18O_depth(count)) < 1e-8
                count2 = count2 + 1;
                d18O_data(count,count2) = d18O_data_temp(n);
            else
                count = count + 1;
                count2 = 1;
                d18O_data(count,count2) = d18O_data_temp(n);
                d18O_depth(count) = depth_d18O_temp(n);
            end
        end
        d18O_depth = d18O_depth(1:count);
        d18O_data = d18O_data(1:count,:);
        d18O_age_info = NaN*ones(count,3);
    else
        d18O_depth = [];
        d18O_data = [];
        d18O_age_info = [];
        MAX = 1;
    end
    
    
    % Gather C14 data:
    path = ['Inputs/',inputFile,'/Records/',data(ll).name,'/C14_data.txt'];
    if strcmp(setting.data_type,'d18O') == 0 && exist(path,'file') == 2
        fileID = fopen(path);
        AA = textscan(fileID,'%s %s %s %s %s %s');
        fclose(fileID);
        C14_data_temp = zeros(length(AA{1})-1,7);
        for k = 1:6
            C14_data_temp(:,k) = str2double(AA{k}(2:end));
        end
        % C14_data_temp = load(path);
        
        C14_depth_temp = C14_data_temp(:,1);
        C14_data_temp(:,2:5) = C14_data_temp(:,2:5);
        C14_data_temp = C14_data_temp(:,2:7);
        
        [~,order] = sort(C14_depth_temp,'ascend');
        C14_depth_temp = C14_depth_temp(order);
        C14_data_temp = C14_data_temp(order,:);
        
        N = length(C14_depth_temp);
        C14_depth = zeros(N,1);
        C14_data = cell(N,1);
        
        count = 1;
        C14_depth(count) = C14_depth_temp(1);
        C14_data{count} = C14_data_temp(1,:);
        for n = 2:N
            if abs(C14_depth(count)-C14_depth_temp(n)) < 1e-8
                C14_data{count} = [C14_data{count};C14_data_temp(n,:)];
            else
                count = count + 1;
                C14_depth(count) = C14_depth_temp(n);
                C14_data{count} = C14_data_temp(n,:);
            end
        end
        C14_depth = C14_depth(1:count);
        C14_data = C14_data(1:count);
    else
        C14_depth = [];
        C14_data = [];
    end
    
    
    % Gather additional depths:
    path = ['Inputs/',inputFile,'/Records/',data(ll).name,'/additional_depths.txt'];
    if exist(path,'file') == 2
        fileID = fopen(path);
        AA = textscan(fileID,'%s');
        fclose(fileID);
        empty_depth = zeros(length(AA{1})-1,1);
        for k = 1:1
            empty_depth(:,k) = str2double(AA{k}(2:end));
        end
    else
        empty_depth = [];
    end
    
    path = ['Inputs/',inputFile,'/Records/',data(ll).name,'/additional_ages.txt'];
    % if strcmp(setting.data_type,'C14') == 0 && exist(path,'file') == 2
    if exist(path,'file') == 2
        fileID = fopen(path);
        BB = textscan(fileID,'%s %s %s %s');
        fclose(fileID);
        AA = zeros(length(BB{1})-1,4);
        for k = 1:4
            AA(:,k) = str2double(BB{k}(2:end));
        end
        
        if size(AA,2) == 4
            add_depth = AA(:,1);
            add_age_info = AA(:,2:4);
            add_data = NaN*ones(length(add_depth),MAX);
            d18O_depth = [d18O_depth;add_depth];
            d18O_data = [d18O_data;add_data];
            d18O_age_info = [d18O_age_info;add_age_info];
            [~,order] = sort(d18O_depth,'ascend');
            d18O_depth = d18O_depth(order);
            d18O_data = d18O_data(order,:);
            d18O_age_info = d18O_age_info(order,:);
            
            N = length(d18O_depth);
            d18O_depth_temp = zeros(N,1);
            d18O_data_temp = zeros(N,MAX);
            d18O_age_info_temp = zeros(N,3);
            d18O_depth_temp(1) = d18O_depth(1);
            d18O_data_temp(1,:) = d18O_data(1,:);
            d18O_age_info_temp(1,:) = d18O_age_info(1,:);
            m = 1;
            for n = 2:N
                if abs(d18O_depth_temp(m)-d18O_depth(n)) < 1e-8
                    if ~isnan(d18O_age_info(n,1))
                        d18O_age_info_temp(m,:) = d18O_age_info(n,:);
                    end
                else
                    m = m + 1;
                    d18O_depth_temp(m) = d18O_depth(n);
                    d18O_data_temp(m,:) = d18O_data(n,:);
                    d18O_age_info_temp(m,:) = d18O_age_info(n,:);
                end
            end
            d18O_depth = d18O_depth_temp(1:m);
            d18O_data = d18O_data_temp(1:m,:);
            d18O_age_info = d18O_age_info_temp(1:m,:);
        end
    end
    
    % Merge them chrononically:
    N = length(d18O_depth) + length(C14_depth) + length(empty_depth);
    merge_depth_temp = [d18O_depth;empty_depth;C14_depth];
    merge_d18O_temp = NaN*ones(N,MAX);
    merge_C14_temp = cell(N,1);
    merge_sugg_age_temp = NaN*ones(N,3);
    
    N1 = length(d18O_depth);
    N2 = length(C14_depth);
    merge_d18O_temp(1:N1,:) = d18O_data;
    merge_sugg_age_temp(1:N1,:) = d18O_age_info;
    merge_C14_temp(end-N2+1:end) = C14_data;
    
    [~,order] = sort(merge_depth_temp,'ascend');
    merge_depth_temp = merge_depth_temp(order);
    merge_d18O_temp = merge_d18O_temp(order,:);
    merge_C14_temp = merge_C14_temp(order);
    merge_sugg_age_temp = merge_sugg_age_temp(order,:);
    
    merge_depth = zeros(N,1);
    merge_d18O = zeros(N,MAX);
    merge_C14 = cell(N,1);
    merge_sugg_age = zeros(N,3);
    
    merge_depth(1) = merge_depth_temp(1);
    merge_d18O(1,:) = merge_d18O_temp(1,:);
    merge_C14{1} = merge_C14_temp{1};
    merge_sugg_age(1,:) = merge_sugg_age_temp(1,:);
    count = 1;
    for n = 2:N
        if abs(merge_depth_temp(n)-merge_depth(count)) < 1e-8
            if ~isnan(merge_d18O_temp(n,1))
                merge_d18O(count,:) = merge_d18O_temp(n,:);
                merge_sugg_age(count,:) = merge_sugg_age_temp(n,:);
            end
            if ~isempty(merge_C14_temp{n})
                merge_C14{count} = merge_C14_temp{n};
            end
        else
            count = count + 1;
            merge_depth(count) = merge_depth_temp(n);
            merge_d18O(count,:) = merge_d18O_temp(n,:);
            merge_C14{count} = merge_C14_temp{n};
            merge_sugg_age(count,:) = merge_sugg_age_temp(n,:);
        end
    end
    merge_depth = merge_depth(1:count);
    merge_d18O = merge_d18O(1:count,:);
    merge_C14 = merge_C14(1:count);
    merge_sugg_age = merge_sugg_age(1:count,:);
    
    data(ll).depth = merge_depth;
    data(ll).d18O = merge_d18O;
    data(ll).radiocarbon = merge_C14;
    data(ll).suggested_age = merge_sugg_age;
    
    path = ['Inputs/',inputFile,'/Records/',data(ll).name,'/initial_ages.txt'];
    if exist(path,'file') == 2
        fileID = fopen(path);
        BB = textscan(fileID,'%s %s');
        fclose(fileID);
        AA = zeros(length(BB{1})-1,2);
        for k = 1:2
            AA(:,k) = str2double(BB{k}(2:end));
        end
        data(ll).initial_age = AA;
    end
end


% record-specific parameters:
for ll = 1:L
    data(ll).phi_I = tran_param(1,:);
    data(ll).phi_C = tran_param(2,:);
    data(ll).phi_M = tran_param(3,:);
    data(ll).phi_E = tran_param(4,:);
    
    path = 'Defaults/setting_core.txt';
    fileID = fopen(path);
    INFO = textscan(fileID,'%s %s');
    fclose(fileID);
    
    data(ll).ACC_MODEL = 'lognormal';
    data(ll).start_depth = str2double(INFO{2}{strcmp(INFO{1},'start_depth:')==1});
    data(ll).end_depth = str2double(INFO{2}{strcmp(INFO{1},'end_depth:')==1});
    data(ll).initial_shift = str2double(INFO{2}{strcmp(INFO{1},'initial_shift:')==1});
    data(ll).initial_scale = str2double(INFO{2}{strcmp(INFO{1},'initial_scale:')==1});
    data(ll).initial_average_sed_rate = str2double(INFO{2}{strcmp(INFO{1},'initial_average_sed_rate:')==1});
    data(ll).islearn_shift = INFO{2}{strcmp(INFO{1},'islearn_shift:')==1};
    data(ll).islearn_scale = INFO{2}{strcmp(INFO{1},'islearn_scale:')==1};
    data(ll).islearn_average_sed_rate = INFO{2}{strcmp(INFO{1},'islearn_average_sed_rate:')==1};
    data(ll).lower_bound = str2double(INFO{2}{strcmp(INFO{1},'lower_bound:')==1});
    data(ll).upper_bound = str2double(INFO{2}{strcmp(INFO{1},'upper_bound:')==1});
    data(ll).min_resolution = str2double(INFO{2}{strcmp(INFO{1},'min_resolution:')==1});
    data(ll).min_resolution_mode = 'absolute';
    data(ll).lower_sedrate = str2double(INFO{2}{strcmp(INFO{1},'lower_sedrate:')==1});
    data(ll).upper_sedrate = str2double(INFO{2}{strcmp(INFO{1},'upper_sedrate:')==1});
    data(ll).is_top_14C_inlier = 'no';
    
    data(ll).SM_BW = str2double(INFO{2}{strcmp(INFO{1},'smoothness_bandwidth:')==1});
    data(ll).PTCL_BW = str2double(INFO{2}{strcmp(INFO{1},'particle_bandwidth:')==1});
    
    path = ['Inputs/',inputFile,'/Records/',data(ll).name,'/setting_core.txt'];
    if exist(path,'file') == 2
        fileID = fopen(path);
        INFO = textscan(fileID,'%s %s');
        fclose(fileID);
        
        if sum(strcmp(INFO{1},'start_depth:')==1) == 1
            data(ll).start_depth = str2double(INFO{2}{strcmp(INFO{1},'start_depth:')==1});
        end
        
        if sum(strcmp(INFO{1},'end_depth:')==1) == 1
            data(ll).end_depth = str2double(INFO{2}{strcmp(INFO{1},'end_depth:')==1});
        end
        
        if sum(strcmp(INFO{1},'initial_shift:')==1) == 1
            data(ll).initial_shift = str2double(INFO{2}{strcmp(INFO{1},'initial_shift:')==1});
        end
        
        if sum(strcmp(INFO{1},'initial_scale:')==1) == 1
            data(ll).initial_scale = str2double(INFO{2}{strcmp(INFO{1},'initial_scale:')==1});
        end
        
        if sum(strcmp(INFO{1},'initial_average_sed_rate:')==1) == 1
            data(ll).initial_average_sed_rate = str2double(INFO{2}{strcmp(INFO{1},'initial_average_sed_rate:')==1});
        end
        
        if sum(strcmp(INFO{1},'islearn_shift:')==1) == 1
            data(ll).islearn_shift = INFO{2}{strcmp(INFO{1},'islearn_shift:')==1};
        end
        
        if sum(strcmp(INFO{1},'islearn_scale:')==1) == 1
            data(ll).islearn_scale = INFO{2}{strcmp(INFO{1},'islearn_scale:')==1};
        end
        
        if sum(strcmp(INFO{1},'islearn_average_sed_rate:')==1) == 1
            data(ll).islearn_average_sed_rate = INFO{2}{strcmp(INFO{1},'islearn_average_sed_rate:')==1};
        end
        
        if sum(strcmp(INFO{1},'lower_bound:')==1) == 1
            data(ll).lower_bound = str2double(INFO{2}{strcmp(INFO{1},'lower_bound:')==1});
        end
        
        if sum(strcmp(INFO{1},'upper_bound:')==1) == 1
            data(ll).upper_bound = str2double(INFO{2}{strcmp(INFO{1},'upper_bound:')==1});
        end
        
        if sum(strcmp(INFO{1},'min_resolution:')==1) == 1
            data(ll).min_resolution = str2double(INFO{2}{strcmp(INFO{1},'min_resolution:')==1});
        end
        
        if sum(strcmp(INFO{1},'lower_sedrate:')==1) == 1
            data(ll).lower_sedrate = str2double(INFO{2}{strcmp(INFO{1},'lower_sedrate:')==1});
        end
        
        if sum(strcmp(INFO{1},'upper_sedrate:')==1) == 1
            data(ll).upper_sedrate = str2double(INFO{2}{strcmp(INFO{1},'upper_sedrate:')==1});
        end
        
        if sum(strcmp(INFO{1},'smoothness_bandwidth:')==1) == 1
            data(ll).SM_BW = str2double(INFO{2}{strcmp(INFO{1},'smoothness_bandwidth:')==1});
        end
        
        if sum(strcmp(INFO{1},'particle_bandwidth:')==1) == 1
            data(ll).PTCL_BW = str2double(INFO{2}{strcmp(INFO{1},'particle_bandwidth:')==1});
        end
    end
    
    path = ['Defaults/Accumulation_Rate_Models/',data(ll).ACC_MODEL,'.txt'];
    data(ll).ACC_MODEL = load(path);
    data(ll).ACC_MODEL(:,2) = interp1(data(ll).ACC_MODEL(:,1),data(ll).ACC_MODEL(:,2),1./data(ll).ACC_MODEL(:,1)).*data(ll).ACC_MODEL(:,1).^(-2);
    ID = (data(ll).ACC_MODEL(:,1)<=1./1.0850);
    CC = data(ll).ACC_MODEL(ID,:);
    DD = [CC(2:end,1)-CC(1:end-1,1),(CC(2:end,2)+CC(1:end-1,2))/2];
    data(ll).ACC_CONTRACTION = sum(DD(:,1).*DD(:,2));
    ID = (data(ll).ACC_MODEL(:,1)<=1./0.9220)&(data(ll).ACC_MODEL(:,1)>1./1.0850);
    CC = data(ll).ACC_MODEL(ID,:);
    DD = [CC(2:end,1)-CC(1:end-1,1),(CC(2:end,2)+CC(1:end-1,2))/2];
    data(ll).ACC_STEADY = sum(DD(:,1).*DD(:,2));
    ID = (data(ll).ACC_MODEL(:,1)>1./0.9220);
    CC = data(ll).ACC_MODEL(ID,:);
    DD = [CC(2:end,1)-CC(1:end-1,1),(CC(2:end,2)+CC(1:end-1,2))/2];
    data(ll).ACC_EXPANSION = sum(DD(:,1).*DD(:,2));
    CC = data(ll).ACC_CONTRACTION + data(ll).ACC_STEADY + data(ll).ACC_EXPANSION;
    data(ll).ACC_CONTRACTION = log(data(ll).ACC_CONTRACTION/CC);
    data(ll).ACC_STEADY = log(data(ll).ACC_STEADY/CC);
    data(ll).ACC_EXPANSION = log(data(ll).ACC_EXPANSION/CC);
    data(ll).ACC_MODEL(:,2) = max(data(ll).ACC_MODEL(:,2),1e-24);
    data(ll).ACC_MODEL(:,2) = log(data(ll).ACC_MODEL(:,2));


    if ~isnan(data(ll).min_resolution)
        data(ll).res = data(ll).min_resolution;
    else
        data(ll).res = inf;
    end
    
    
    if isnan(data(ll).upper_sedrate)
        data(ll).upper_sedrate = inf;
    end

    if isnan(data(ll).lower_sedrate) || (data(ll).lower_sedrate<0)
        data(ll).lower_sedrate = 0;
    end
    
    if isnan(data(ll).initial_scale)
        data(ll).scale = 1;
    else
        data(ll).scale = data(ll).initial_scale;
    end
    
    if isnan(data(ll).initial_shift)
        index = (~isnan(data(ll).d18O));
        AA = data(ll).d18O(index);
        BB = target.stack(:,2)*data(ll).scale;
        QQ = max(AA) + min(AA);
        PP = max(BB) + min(BB);
        data(ll).shift = (QQ-PP)/2;
        if isempty(data(ll).shift)
            data(ll).shift = 0;
        end
    else
        data(ll).shift = data(ll).initial_shift;
    end
    
    if ~isnan(data(ll).lower_bound)
        if strcmp(setting.data_type,'C14')
            data(ll).min = max(data(ll).lower_bound,0);
        else
            data(ll).min = max(data(ll).lower_bound,target.stack(1,1));
        end
    else
        if strcmp(setting.data_type,'C14')
            data(ll).min = 0;
        else
            data(ll).min = target.stack(1,1);
        end
    end
    
    if ~isnan(data(ll).upper_bound)
        if strcmp(setting.data_type,'C14')
            % data(ll).max = min(data(ll).upper_bound,60);
            data(ll).max = data(ll).upper_bound;
        else
            data(ll).max = min(data(ll).upper_bound,target.stack(end,1));
        end
    else
        if strcmp(setting.data_type,'C14')
            data(ll).max = 60;
        else
            data(ll).max = target.stack(end,1);
        end
    end
    [QR1,QR2] = max(data(ll).suggested_age(:,1));
    if ~isnan(QR1)
        data(ll).max = max(data(ll).max,data(ll).suggested_age(QR2,1)+3*data(ll).suggested_age(QR2,2));
    end
    
    if isnan(data(ll).initial_average_sed_rate)
        T = size(target.stack,1);
        N = length(data(ll).depth);
        age_info = data(ll).suggested_age;
        depth = data(ll).depth;
        C14 = data(ll).radiocarbon;
        
        st = [inf;0];
        ed = [-inf;0];
        for n = 1:N
            if ~isnan(age_info(n,1)) || ~isempty(C14{n})
                st(1) = min(st(1),n);
                ed(1) = max(ed(1),n);
            end
        end
        
        if ~isinf(st(1))
            if ~isnan(age_info(st(1),1))
                st(2) = 1;
            end
        end
        
        if ~isinf(ed(1))
            if ~isnan(age_info(ed(1),1))
                ed(2) = 1;
            end
        end
        
        st_age = target.stack(1,1);
        ed_age = target.stack(end,1);
        if ~isinf(ed(1))
            if ed(2) == 1
                ed_age = age_info(ed(1),1);
            else
                [AA,order] = unique(target.cal_curve{C14{ed(1)}(1,5)}(:,2));
                BB = target.cal_curve{C14{ed(1)}(1,5)}(order,1);
                ed_age = interp1(AA,BB,min(max(C14{ed(1)}(1,1)-C14{ed(1)}(1,3),target.cal_curve{C14{ed(1)}(1,5)}(1,2)),target.cal_curve{C14{ed(1)}(1,5)}(end,2)));
            end
            if ~isinf(st(1)) && st(1) < ed(1)
                if st(2) == 1
                    st_age = age_info(st(1),1);
                else
                    [AA,order] = unique(target.cal_curve{C14{st(1)}(1,5)}(:,2));
                    BB = target.cal_curve{C14{st(1)}(1,5)}(order,1);
                    st_age = interp1(AA,BB,min(max(C14{st(1)}(1,1)-C14{st(1)}(1,3),target.cal_curve{C14{st(1)}(1,5)}(1,2)),target.cal_curve{C14{st(1)}(1,5)}(end,2)));
                end
            else
                st(1) = 1;
            end
        else
            st(1) = 1;
            ed(1) = N;
        end
        if st(1) == ed(1)
            data(ll).R = [target.stack(:,1),(target.stack(end,1)-target.stack(1,1))/(depth(end)-depth(1))*ones(T,1)];
        else
            if strcmp(setting.data_type,'C14')
                RT = (data(ll).min:0.5:data(ll).max)';
                if RT(end) < data(ll).max
                    RT = [RT;data(ll).max];
                end
                data(ll).R = [RT,(ed_age-st_age)/(depth(ed(1))-depth(st(1)))*ones(size(RT,1),1)];
            else
                data(ll).R = [target.stack(:,1),(ed_age-st_age)/(depth(ed(1))-depth(st(1)))*ones(T,1)];
            end
        end
    else
        data(ll).R = [target.stack(:,1),1/initial_average_sed_rate*ones(T,1)];
    end
    
    if strcmp(data(ll).is_top_14C_inlier,'yes')
        N = length(data(ll).radiocarbon);
        
        DET = 1;
        for n = 1:N
            if DET == 1 && ~isempty(data(ll).radiocarbon{n})
                DET = 0;
                data(ll).radiocarbon{n}(:,end) = 1;
            end
        end
    end
    
    if ~isnan(data(ll).start_depth)
        ID = (data(ll).depth>=data(ll).start_depth);
        data(ll).depth = data(ll).depth(ID,:);
        data(ll).d18O = data(ll).d18O(ID,:);
        data(ll).radiocarbon = data(ll).radiocarbon(ID,:);
        data(ll).suggested_age = data(ll).suggested_age(ID,:);
    end
    
    if ~isnan(data(ll).end_depth)
        ID = (data(ll).depth<=data(ll).end_depth);
        data(ll).depth = data(ll).depth(ID,:);
        data(ll).d18O = data(ll).d18O(ID,:);
        data(ll).radiocarbon = data(ll).radiocarbon(ID,:);
        data(ll).suggested_age = data(ll).suggested_age(ID,:);
    end
end

data = rmfield(data,'initial_shift');
data = rmfield(data,'initial_scale');
data = rmfield(data,'initial_average_sed_rate');
data = rmfield(data,'lower_bound');
data = rmfield(data,'upper_bound');


% add depths:
for ll = 1:L
    depth = data(ll).depth;
    d18O = data(ll).d18O;
    C14 = data(ll).radiocarbon;
    add_age = data(ll).suggested_age;
    
    N = length(depth);
    
    if strcmp(data(ll).min_resolution_mode,'off') == 1
        d = inf;
    elseif strcmp(data(ll).min_resolution_mode,'absolute') == 1
        d = data(ll).res; % in meters
    elseif strcmp(data(ll).min_resolution_mode,'relative') == 1
        D = depth(end)-depth(1);
        d = D*data(ll).res; % proporational to the core length
    end
    
    aug_depth = [];
    m = 0;
    for n = 1:N-1
        r = floor(abs(depth(n)-depth(n+1))/d);
        if r > 0
            A = linspace(depth(n),depth(n+1),r+2);
            aug_depth(m+1:m+r) = A(2:end-1);
            m = m + r;
        end
    end
    aug_depth = aug_depth';
    
    M = length(aug_depth);
    
    depth = [depth;aug_depth];
    d18O = [d18O;NaN*ones(M,size(d18O,2))];
    C14 = [C14;cell(M,1)];
    add_age = [add_age;NaN*ones(M,3)];
    
    [~,order] = sort(depth,'ascend');
    depth = depth(order,:);
    d18O = d18O(order,:);
    C14 = C14(order,:);
    add_age = add_age(order,:);
    
    data(ll).depth = depth;
    data(ll).d18O = d18O;
    data(ll).radiocarbon = C14;
    data(ll).suggested_age = add_age;
end


for ll = 1:L
    data(ll).suggested_age(:,3) = 1 - data(ll).suggested_age(:,3);
end



% truncate records:
data_ps = data;

for ll = 1:L
    N = length(data(ll).depth);
    
    dT = 1/mean(data(ll).R(:,2));
    
    depth = zeros(N,1);
    d18O = zeros(size(data(ll).d18O));
    radiocarbon = cell(N,1);
    sugg_age = zeros(N,3);
    index = zeros(N,1);
    
    depth(1) = data(ll).depth(1);
    d18O(1,:) = data(ll).d18O(1,:);
    radiocarbon{1} = data(ll).radiocarbon{1};
    sugg_age(1,:) = data(ll).suggested_age(1,:);
    index(1) = 1;
    d = 0;
    m = 1;
    for n = 2:N-1
        d = d + (data(ll).depth(n)-data(ll).depth(n-1));
        if ~isempty(data(ll).radiocarbon{n}) || ~isnan(data(ll).suggested_age(n,1))
            m = m + 1;
            depth(m) = data(ll).depth(n);
            d18O(m,:) = data(ll).d18O(n,:);
            radiocarbon{m} = data(ll).radiocarbon{n};
            sugg_age(m,:) = data(ll).suggested_age(n,:);
            index(m) = n;
            d = 0;
        else
            if ~isempty(data(ll).initial_age) && sum(abs(data(ll).initial_age(:,1)-data(ll).depth(n))<1e-6)>0
                m = m + 1;
                depth(m) = data(ll).depth(n);
                d18O(m,:) = data(ll).d18O(n,:);
                radiocarbon{m} = data(ll).radiocarbon{n};
                sugg_age(m,:) = data(ll).suggested_age(n,:);
                index(m) = n;
                d = 0;
            else
                if d > dT && ~isnan(data(ll).d18O(n,1))
                    m = m + 1;
                    depth(m) = data(ll).depth(n);
                    d18O(m,:) = data(ll).d18O(n,:);
                    radiocarbon{m} = data(ll).radiocarbon{n};
                    sugg_age(m,:) = data(ll).suggested_age(n,:);
                    index(m) = n;
                    d = 0;
                end
            end
        end
    end
    m = m + 1;
    
    depth(m) = data(ll).depth(N);
    d18O(m,:) = data(ll).d18O(N,:);
    radiocarbon{m} = data(ll).radiocarbon{N};
    sugg_age(m,:) = data(ll).suggested_age(N,:);
    index(m) = N;
    
    data_ps(ll).depth = depth(1:m);
    data_ps(ll).d18O = d18O(1:m,:);
    data_ps(ll).radiocarbon = radiocarbon(1:m);
    data_ps(ll).suggested_age = sugg_age(1:m,:);
    data_ps(ll).index = index(1:m);
end


end