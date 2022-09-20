function AgeVsDepth(results_path)

results = load([results_path,'/results.mat']);

data = results.summary;
setting = results.setting_alignment;

data_type = setting.data_type;
if strcmp(data_type,'both') == 1
    data_type = 'dual proxy';
end

path = [results_path,'/figures/age_vs_depth'];
if exist(path,'dir') ~= 7
    mkdir(path);
end

if strcmp(setting.data_type,'both') || strcmp(setting.data_type,'C14')
    CI = results.CI_C14;
end

clear results;

M = length(data);

for m = 1:M
    Age = data(m).age_samples;
    Depth = data(m).depth;
    
    N = length(Depth);
    LMU = zeros(N,3);
    
    for n = 1:N
        LMU(n,1) = quantile(Age(n,:),0.025);
        LMU(n,2) = quantile(Age(n,:),0.5);
        LMU(n,3) = quantile(Age(n,:),0.975);
    end
    
    MAX = LMU(N,3);
    MAX = ceil(MAX/10)*10;
    
    TYPE = zeros(N,1);
    for n = 1:N
        if isnan(data(m).d18O(n))
            if isempty(data(m).radiocarbon{n})
                % additional depths have label 1
                TYPE(n) = 1;
            else
                % radiocarbon depths have label 2
                TYPE(n) = 2;
            end
        else
            if isempty(data(m).radiocarbon{n})
                % d18O depths have label 3
                TYPE(n) = 3;
            else
                % depths of both kindes of data have label 4
                TYPE(n) = 4;
            end
        end
    end
    
    h = zeros(1,1);
    fig = figure;
    
    subplot(2,3,[1,2]);
    hold on;
    tt = data(m).name;
    tt(tt=='_') = '-';
    title(tt,'FontSize',16);
    xx = [LMU(:,1);flipud(LMU(:,3))];
    yy = [data(m).depth;flipud(data(m).depth)];
    patch(xx,yy,1,'FaceColor','k','FaceAlpha',0.1,'EdgeColor','none');
    h(1) = plot(LMU(:,2),data(m).depth,'--','LineWidth',2,'Color','k');
    if strcmp(setting.data_type,'both') || strcmp(setting.data_type,'C14')
        if isempty(CI(m).depth) == 0
            q = 2;
            for j = 1:length(CI(m).depth)
                h(2) = plot([CI(m).lower_95(j),CI(m).upper_95(j)],[CI(m).depth(j),CI(m).depth(j)],'b','LineWidth',2);
            end
        else
            q = 1;
        end
    else
        q = 1;
    end
    
    DET = 0;
    for n = 1:N
        if ~isnan(data(m).additional_ages(n,1))
            DET = 1;
            mu = data(m).additional_ages(n,1);
            sig = data(m).additional_ages(n,2);
            if data(m).additional_ages(n,3) == 0
                h(q+1) = plot([max(mu-1.96*sig,0),min(mu+1.96*sig,MAX)],[data(m).depth(n),data(m).depth(n)],'--m','LineWidth',2);
            elseif data(m).additional_ages(n,3) == 1
                h(q+1) = plot([max(mu-sig,0),min(mu+sig,MAX)],[data(m).depth(n),data(m).depth(n)],'m','LineWidth',2);
            end
        end
    end
    
    if DET == 1
        if q == 1
            list = {data_type,'Age Guess'};
        elseif q == 2
            list = {data_type,'Individual 14C','Age Guess'};
        end
        q = q + 1;
    else
        if q == 1
            list = {data_type};
        elseif q == 2
            list = {data_type,'Individual 14C'};
        end
    end
    
    if sum(TYPE == 3) > 0
        q = q + 1;
        h(q) = plot(zeros(sum(TYPE == 3),1),data(m).depth(TYPE == 3,:),'>r');
        list{q} = 'd18O';
    end
    if sum(TYPE == 2) > 0
        q = q + 1;
        h(q) = plot(zeros(sum(TYPE == 2),1),data(m).depth(TYPE == 2,:),'>c');
        list{q} = '14C';
    end
    if sum(TYPE == 4) > 0
        q = q + 1;
        h(q) = plot(zeros(sum(TYPE == 4),1),data(m).depth(TYPE == 4,:),'>g');
        list{q} = 'Both';
    end
    if sum(TYPE == 1) > 0
        q = q + 1;
        h(q) = plot(zeros(sum(TYPE == 1),1),data(m).depth(TYPE == 1,:),'>k');
        list{q} = 'None';
    end
    
    xlabel('age (ky)','FontSize',12);
    ylabel('depth (m)','FontSize',12);
    xlim([0 MAX]);
    ylim([0 data(m).depth(end)*1.05]);
    legend(h,list,'Location','SouthEast','FontSize',12);
    
    
    
    subplot(2,3,[4,5]);
    hold on;
    for t = 1:size(Age,2)
        h(1) = plot(Age(:,t),data(m).depth,'.k');
    end
    if strcmp(setting.data_type,'both') || strcmp(setting.data_type,'C14')
        if isempty(CI(m).depth) == 0
            q = 2;
            for j = 1:length(CI(m).depth)
                h(2) = plot([CI(m).lower_95(j),CI(m).upper_95(j)],[CI(m).depth(j),CI(m).depth(j)],'c','LineWidth',2);
            end
        else
            q = 1;
        end
    else
        q = 1;
    end
    
    DET = 0;
    for n = 1:N
        if ~isnan(data(m).additional_ages(n,1))
            DET = 1;
            mu = data(m).additional_ages(n,1);
            sig = data(m).additional_ages(n,2);
            if data(m).additional_ages(n,3) == 0
                h(q+1) = plot([max(mu-1.96*sig,0),min(mu+1.96*sig,MAX)],[data(m).depth(n),data(m).depth(n)],'--m','LineWidth',2);
            elseif data(m).additional_ages(n,3) == 1
                h(q+1) = plot([max(mu-sig,0),min(mu+sig,MAX)],[data(m).depth(n),data(m).depth(n)],'m','LineWidth',2);
            end
        end
    end
    
    if DET == 1
        if q == 1
            list = {data_type,'Age Guess'};
        elseif q == 2
            list = {data_type,'Individual 14C','Age Guess'};
        end
        q = q + 1;
    else
        if q == 1
            list = {data_type};
        elseif q == 2
            list = {data_type,'Individual 14C'};
        end
    end
    
    if sum(TYPE == 3) > 0
        q = q + 1;
        h(q) = plot(zeros(sum(TYPE == 3),1),data(m).depth(TYPE == 3,:),'>r');
        list{q} = 'd18O';
    end
    if sum(TYPE == 2) > 0
        q = q + 1;
        h(q) = plot(zeros(sum(TYPE == 2),1),data(m).depth(TYPE == 2,:),'>c');
        list{q} = '14C';
    end
    if sum(TYPE == 4) > 0
        q = q + 1;
        h(q) = plot(zeros(sum(TYPE == 4),1),data(m).depth(TYPE == 4,:),'>g');
        list{q} = 'Both';
    end
    if sum(TYPE == 1) > 0
        q = q + 1;
        h(q) = plot(zeros(sum(TYPE == 1),1),data(m).depth(TYPE == 1,:),'>k');
        list{q} = 'None';
    end
    
    xlabel('age (ky)','FontSize',12);
    ylabel('depth (m)','FontSize',12);
    xlim([0 MAX]);
    ylim([0 data(m).depth(end)*1.05]);
    legend(h,list,'Location','SouthEast','FontSize',12);
    
    
    QQ = LMU(:,2);
    
    LMU(:,1) = LMU(:,1) - QQ;
    LMU(:,2) = LMU(:,2) - QQ;
    LMU(:,3) = LMU(:,3) - QQ;
    
    MAX = 6;
    MIN = -6;
    
    subplot(2,3,3);
    hold on;
    xx = [LMU(:,1);flipud(LMU(:,3))];
    yy = [data(m).depth;flipud(data(m).depth)];
    patch(xx,yy,1,'FaceColor','k','FaceAlpha',0.1,'EdgeColor','none');
    plot(LMU(:,2),data(m).depth,'--','LineWidth',2,'Color','k');
    xlabel('relative age (ky)','FontSize',12);
    xlim([MIN MAX]);
    ylim([0 data(m).depth(end)*1.05]);
    
    set(fig,'Position',[10 10 1200 1100]);
    movegui(fig,'center');
    
    path = [results_path,'/figures/age_vs_depth/',data(m).name,'.fig'];
    savefig(fig,path);
end

close all;


end