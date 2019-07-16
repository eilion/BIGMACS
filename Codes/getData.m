function [data] = getData(record_path,data_type)

core_title = textread(record_path,'%s');

L = length(core_title);

data = struct('name',cell(L,1),'depth',cell(L,1),'d18O',cell(L,1),'radiocarbon',cell(L,1),'suggested_age',cell(L,1));

for ll = 1:L
    data(ll).name = core_title{ll};
    
    % Gather d18O data:
    if strcmp(data_type,'C14') == 0
        path = ['Cores/',data(ll).name,'/d18O_data.txt'];
        d18O_data_temp = load(path);
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
    if strcmp(data_type,'d18O') == 0
        path = ['Cores/',data(ll).name,'/C14_data.txt'];
        C14_data_temp = load(path);
        C14_depth_temp = C14_data_temp(:,1)/100;
        C14_data_temp(:,2:5) = C14_data_temp(:,2:5)/1000;
        C14_data_temp = C14_data_temp(:,2:6);
        
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
    path = ['Cores/',data(ll).name,'/additional_ages.txt'];
    if strcmp(data_type,'C14') == 0 && exist(path,'file') == 2
        AA = load(path);
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
                    if isnan(d18O_age_info(n,1)) == 0
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
    N = length(d18O_depth) + length(C14_depth);
    merge_depth_temp = [d18O_depth;C14_depth];
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
            if isnan(merge_d18O_temp(n,1)) == 0
                merge_d18O(count,:) = merge_d18O_temp(n,:);
                merge_sugg_age(count,:) = merge_sugg_age_temp(n,:);
            end
            if isempty(merge_C14_temp{n}) == 0
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
end


end

