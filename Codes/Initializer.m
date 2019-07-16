function [core_param,global_param] = Initializer(data,stack,cal_curve,inputFile_path,data_type)

L = length(data);
T = size(stack,1);

path = [inputFile_path,'/transition_parameter.txt'];
if exist(path,'file') ~= 2
    path = 'Defaults/transition_parameter.txt';
end
tran_param = load(path);

global_param = getHyperparameters(inputFile_path);
global_param.tran_param = tran_param;

core_param = struct('name',cell(L,1),'scale',cell(L,1),'shift',cell(L,1),'R',cell(L,1),'phi_I',cell(L,1),'phi_C',cell(L,1),'phi_M',cell(L,1),'phi_E',cell(L,1),'min',cell(L,1),'max',cell(L,1),'res',cell(L,1),'lower_sedrate',cell(L,1),'upper_sedrate',cell(L,1));

for ll = 1:L
    core_param(ll).name = data(ll).name;
    
    core_param(ll).phi_I = tran_param(1,:);
    core_param(ll).phi_C = tran_param(2,:);
    core_param(ll).phi_M = tran_param(3,:);
    core_param(ll).phi_E = tran_param(4,:);
    
    [initial_shift,initial_scale,initial_average_sed_rate,~,~,~,lower_bound,upper_bound,min_resolution,lower_sedrate,upper_sedrate] = getSetting_Core(data(ll).name);
    
    if isnan(min_resolution) == 0
        core_param(ll).res = min_resolution;
    else
        core_param(ll).res = inf;
    end
    
    if isinf(upper_sedrate) == 1 || isnan(upper_sedrate) == 1
        core_param(ll).lower_sedrate = -inf;
    else
        core_param(ll).lower_sedrate = 1/upper_sedrate;
    end
    if isinf(lower_sedrate) == 1 || isnan(lower_sedrate) == 1
        core_param(ll).upper_sedrate = inf;
    else
        core_param(ll).upper_sedrate = 1/lower_sedrate;
    end
    
    if isnan(initial_scale) == 1
        core_param(ll).scale = 1;
    else
        core_param(ll).scale = initial_scale;
    end
    
    if isnan(initial_shift) == 1
        index = (isnan(data(ll).d18O) == 0);
        AA = data(ll).d18O(index);
        BB = stack(:,2)*core_param(ll).scale;
        QQ = max(AA) + min(AA);
        PP = max(BB) + min(BB);
        core_param(ll).shift = (QQ-PP)/2;
        if isempty(core_param(ll).shift) == 1
            core_param(ll).shift = 0;
        end
    else
        core_param(ll).shift = initial_shift;
    end
    
    if isnan(lower_bound) == 0
        if strcmp(data_type,'C14') == 1
            core_param(ll).min = max(lower_bound,0);
        else
            core_param(ll).min = max(lower_bound,stack(1,1));
        end
    else
        if strcmp(data_type,'C14') == 1
            core_param(ll).min = 0;
        else
            core_param(ll).min = stack(1,1);
        end
    end
    
    if isnan(upper_bound) == 0
        if strcmp(data_type,'C14') == 1
            core_param(ll).max = min(upper_bound,60);
        else
            core_param(ll).max = min(upper_bound,stack(end,1));
        end
    else
        if strcmp(data_type,'C14') == 1
            core_param(ll).max = 60;
        else
            core_param(ll).max = stack(end,1);
        end
    end
    
    if isnan(initial_average_sed_rate) == 1
        N = length(data(ll).depth);
        age_info = data(ll).suggested_age;
        depth = data(ll).depth;
        C14 = data(ll).radiocarbon;
        
        st = [inf;0];
        ed = [-inf;0];
        for n = 1:N
            if isnan(age_info(n,1)) == 0 || isempty(C14{n}) == 0
                st(1) = min(st(1),n);
                ed(1) = max(ed(1),n);
            end
        end
        
        if isinf(st(1)) == 0
            if isnan(age_info(st(1),1)) == 0
                st(2) = 1;
            end
        end
        
        if isinf(ed(1)) == 0
            if isnan(age_info(ed(1),1)) == 0
                ed(2) = 1;
            end
        end
        
        st_age = stack(1,1);
        ed_age = stack(end,1);
        if isinf(ed(1)) == 0
            if ed(2) == 1
                ed_age = age_info(ed(1),1);
            else
                ed_age = interp1(cal_curve{C14{ed(1)}(1,5)}(:,2),cal_curve{C14{ed(1)}(1,5)}(:,1),min(max(C14{ed(1)}(1,1)-C14{ed(1)}(1,3),cal_curve{C14{ed(1)}(1,5)}(1,2)),cal_curve{C14{ed(1)}(1,5)}(end,2)));
            end
            if isinf(st(1)) == 0 && st(1) < ed(1)
                if st(2) == 1
                    st_age = age_info(st(1),1);
                else
                    st_age = interp1(cal_curve{C14{st(1)}(1,5)}(:,2),cal_curve{C14{st(1)}(1,5)}(:,1),min(max(C14{st(1)}(1,1)-C14{st(1)}(1,3),cal_curve{C14{st(1)}(1,5)}(1,2)),cal_curve{C14{st(1)}(1,5)}(end,2)));
                end
            else
                st(1) = 1;
            end
        else
            st(1) = 1;
            ed(1) = N;
        end
        core_param(ll).R = [stack(:,1),(ed_age-st_age)/(depth(ed(1))-depth(st(1)))*ones(T,1)];
    else
        core_param(ll).R = [stack(:,1),1/initial_average_sed_rate*ones(T,1)];
    end
end


end