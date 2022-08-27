function AgeVsSedRate(results_path)

results = load([results_path,'/results.mat']);

path = [results_path,'/figures/age_vs_sed_rate'];
if exist(path,'dir') ~= 7
    mkdir(path);
end

samples = results.samples;

L = 1;

clear results;

cc = zeros(1,3);

M = length(samples);

SAMPLES = cell(L,1);
data = cell(L,1);
data_type = cell(L,1);

results = load([results_path,'/results.mat']);

SAMPLES{1} = results.samples;
data{1} = results.data;
data_type{1} = results.setting.data_type;
if strcmp(data_type{1},'both') == 1
    data_type{1} = 'dual proxy';
end
clear results;

for m = 1:M
    age = cell(L,1);
    log_sed_rate = cell(L,1);
    
    MEANs = cell(L,1);
    As = cell(L,1);
    Bs = cell(L,1);
    Cx = cell(L,1);
    Cy = cell(L,1);
    Angle = cell(L,1);
    
    MAX = -inf;
    
    for ll = 1:L
        Age = SAMPLES{ll}(m).ages;
        Depth = SAMPLES{ll}(m).depth;
        
        % age{ll} = (Age(2:end,:)+Age(1:end-1,:))/2;
        age{ll} = Age(2:end,:);
        log_sed_rate{ll} = log((Depth(2:end)-Depth(1:end-1))./(Age(2:end,:)-Age(1:end-1,:)));
        
        MAX = max([MAX,quantile(Age(end,:),0.975)]);
        
        N = size(age{ll},1);
        
        MEANs{ll} = zeros(2,N);
        As{ll} = zeros(N,1);
        Bs{ll} = zeros(N,1);
        Cx{ll} = zeros(N,1);
        Cy{ll} = zeros(N,1);
        Angle{ll} = zeros(N,1);
        
        for n = 1:N
            data_pair = [age{ll}(n,:);log_sed_rate{ll}(n,:)];
            MEAN = sum(data_pair,2)/size(data_pair,2);
            MEANs{ll}(:,n) = MEAN;
            COV = (data_pair-MEAN)*(data_pair-MEAN)'/size(data_pair,2);
            [V,D] = eig(COV);
            if D(1,1) < D(2,2)
                VV = V(:,2);
                DD = [D(2,2);D(1,1)];
            else
                VV = V(:,1);
                DD = [D(1,1);D(2,2)];
            end
            As{ll}(n) = sqrt(5.991*DD(1));
            Bs{ll}(n) = sqrt(5.991*DD(2));
            Cx{ll}(n) = MEAN(1);
            Cy{ll}(n) = MEAN(2);
            Angle{ll}(n) = atan(VV(2)/VV(1));
        end
    end
    MAX = ceil(MAX/10)*10;
    
    h = zeros(L,1);
    fig = figure;
    subplot(2,1,1);
    hold on;
    tt = SAMPLES{1}(m).name;
    tt(tt=='_') = '-';
    title(tt,'FontSize',16);
    for ll = 1:L
        plot(data{ll}(m).R(:,1),1./data{ll}(m).R(:,2),'Linewidth',2,'Color',cc(ll,:));
        Elliptic_Confidence_Region(As{ll},Bs{ll},Cx{ll},Cy{ll},Angle{ll},cc(ll,:),'sed_rate');
    end
    for ll = 1:L
        h(ll) = plot(MEANs{ll}(1,:),exp(MEANs{ll}(2,:)),':*','Linewidth',2,'Color',cc(ll,:));
    end
    xlim([0 MAX]);
    xlabel('age (ky)','FontSize',12);
    ylabel('sedimentation rate (m/ky)','FontSize',12);
    legend(h,data_type,'Location','EastOutside');
    
    age = cell(L,1);
    log_sed_rate = cell(L,1);
    
    MEANs = cell(L,1);
    As = cell(L,1);
    Bs = cell(L,1);
    Cx = cell(L,1);
    Cy = cell(L,1);
    Angle = cell(L,1);
    
    MAX = -inf;
    
    for ll = 1:L
        Age = SAMPLES{ll}(m).ages;
        Depth = SAMPLES{ll}(m).depth;
        R = data{ll}(m).R;
        
        % age{ll} = (Age(2:end,:)+Age(1:end-1,:))/2;
        age{ll} = Age(2:end,:);
        RR = interp1(R(:,1),R(:,2),age{ll});
        log_sed_rate{ll} = log((Depth(2:end)-Depth(1:end-1))./(Age(2:end,:)-Age(1:end-1,:))) + log(RR);
        
        MAX = max([MAX,quantile(Age(end,:),0.975)]);
        
        N = size(age{ll},1);
        
        MEANs{ll} = zeros(2,N);
        As{ll} = zeros(N,1);
        Bs{ll} = zeros(N,1);
        Cx{ll} = zeros(N,1);
        Cy{ll} = zeros(N,1);
        Angle{ll} = zeros(N,1);
        
        for n = 1:N
            data_pair = [age{ll}(n,:);log_sed_rate{ll}(n,:)];
            MEAN = sum(data_pair,2)/size(data_pair,2);
            MEANs{ll}(:,n) = MEAN;
            COV = (data_pair-MEAN)*(data_pair-MEAN)'/size(data_pair,2);
            [V,D] = eig(COV);
            if D(1,1) < D(2,2)
                VV = V(:,2);
                DD = [D(2,2);D(1,1)];
            else
                VV = V(:,1);
                DD = [D(1,1);D(2,2)];
            end
            As{ll}(n) = sqrt(5.991*DD(1));
            Bs{ll}(n) = sqrt(5.991*DD(2));
            Cx{ll}(n) = MEAN(1);
            Cy{ll}(n) = MEAN(2);
            Angle{ll}(n) = atan(VV(2)/VV(1));
        end
    end
    MAX = ceil(MAX/10)*10;
    
    h = zeros(L,1);
    subplot(2,1,2);
    hold on;
    tt = SAMPLES{1}(m).name;
    tt(tt=='_') = '-';
    title(tt,'FontSize',16);
    for ll = 1:L
        Elliptic_Confidence_Region(As{ll},Bs{ll},Cx{ll},Cy{ll},Angle{ll},cc(ll,:),'log_sed_rate');
    end
    for ll = 1:L
        h(ll) = plot(MEANs{ll}(1,:),MEANs{ll}(2,:),':*','Linewidth',2,'Color',cc(ll,:));
    end
    for ll = 1:L
        if isfield(data{ll},'lower_sedrate') == 1
            if ~isnan(data{ll}(m).lower_sedrate) && ~isinf(data{ll}(m).lower_sedrate)
                plot([0,MAX],[log(data{ll}(m).lower_sedrate),log(data{ll}(m).lower_sedrate)],'--','Linewidth',2,'Color',cc(ll,:));
            end
        end
        if isfield(data{ll},'upper_sedrate') == 1
            if ~isnan(data{ll}(m).upper_sedrate) && ~isinf(data{ll}(m).upper_sedrate)
                plot([0,MAX],[log(data{ll}(m).upper_sedrate),log(data{ll}(m).upper_sedrate)],'--','Linewidth',2,'Color',cc(ll,:));
            end
        end
    end
    xlim([0 MAX]);
    xlabel('age (ky)','FontSize',12);
    ylabel('normalized log-sedimentation rate','FontSize',12);
    legend(h,data_type,'Location','EastOutside');
    
    set(fig,'Position',[10 10 1200 1100]);
    movegui(fig,'center');
    
    path = [results_path,'/figures/age_vs_sed_rate/',SAMPLES{1}(m).name,'.fig'];
    savefig(fig,path);
end

close all;


end