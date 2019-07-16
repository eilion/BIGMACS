function AgeVsDepth(results_path)

L = length(results_path);

if L == 1
    results = load([results_path{1},'/results.mat']);
    results = results.results;
    data = results.data;
    samples = results.samples;
    data_type = results.type;
    if strcmp(data_type,'both') == 1
        data_type = 'dual proxy';
    end
    
    CI = getCI_C14(data);
    
    M = length(samples);
    
    for m = 1:M
        Age = samples(m).ages;
        Depth = samples(m).depth;
        
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
            if isnan(data(m).d18O(n)) == 1
                if isempty(data(m).radiocarbon{n}) == 1
                    % additional depths have label 1
                    TYPE(n) = 1;
                else
                    % radiocarbon depths have label 2
                    TYPE(n) = 2;
                end
            else
                if isempty(data(m).radiocarbon{n}) == 1
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
        if isempty(CI(m).depth) == 0
            q = 2;
            for j = 1:length(CI(m).depth)
                h(2) = plot([CI(m).lower(j),CI(m).upper(j)],[CI(m).depth(j),CI(m).depth(j)],'c','LineWidth',2);
            end
        else
            q = 1;
        end
        
        DET = 0;
        for n = 1:N
            if isnan(data(m).suggested_age(n,1)) == 0
                DET = 1;
                mu = data(m).suggested_age(n,1);
                sig = data(m).suggested_age(n,2);
                if data(m).suggested_age(n,3) == 0
                    h(q+1) = plot([max(mu-1.96*sig,0),min(mu+1.96*sig,MAX)],[data(m).depth(n),data(m).depth(n)],'--m','LineWidth',2);
                elseif data(m).suggested_age(n,3) == 1
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
        if isempty(CI(m).depth) == 0
            q = 2;
            for j = 1:length(CI(m).depth)
                h(2) = plot([CI(m).lower(j),CI(m).upper(j)],[CI(m).depth(j),CI(m).depth(j)],'c','LineWidth',2);
            end
        else
            q = 1;
        end
        
        DET = 0;
        for n = 1:N
            if isnan(data(m).suggested_age(n,1)) == 0
                DET = 1;
                mu = data(m).suggested_age(n,1);
                sig = data(m).suggested_age(n,2);
                if data(m).suggested_age(n,3) == 0
                    h(q+1) = plot([max(mu-1.96*sig,0),min(mu+1.96*sig,MAX)],[data(m).depth(n),data(m).depth(n)],'--m','LineWidth',2);
                elseif data(m).suggested_age(n,3) == 1
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
    end
elseif L > 1
    results = load([results_path{1},'/results.mat']);
    results = results.results;
    samples = results.samples;
    
    COLOR = jet(101);
    ss = (1:101)';
    tt = linspace(1,101,L)';
    cc = zeros(L,3);
    cc(:,1) = interp1(ss,COLOR(:,1),tt);
    cc(:,2) = interp1(ss,COLOR(:,2),tt);
    cc(:,3) = interp1(ss,COLOR(:,3),tt);
    
    M = length(samples);
    
    DATA = cell(L,1);
    SAMPLES = cell(L,1);
    data_type = cell(L,1);
    CI = cell(L,1);
    
    for ll = 1:L
        results = load([results_path{ll},'/results.mat']);
        results = results.results;
        DATA{ll} = results.data;
        SAMPLES{ll} = results.samples;
        data_type{ll} = results.type;
        if strcmp(data_type{ll},'both') == 1
            data_type{ll} = 'dual proxy';
        end
        CI{ll} = getCI_C14(DATA{ll});
    end
    
    for m = 1:M
        LMU = cell(L,1);
        MAX = zeros(L,1);
        DEPTH = cell(L,1);
        
        for ll = 1:L
            Age = SAMPLES{ll}(m).ages;
            Depth = SAMPLES{ll}(m).depth;
            
            N = length(Depth);
            LMU{ll} = zeros(N,3);
            
            for n = 1:N
                LMU{ll}(n,1) = quantile(Age(n,:),0.025);
                LMU{ll}(n,2) = quantile(Age(n,:),0.5);
                LMU{ll}(n,3) = quantile(Age(n,:),0.975);
            end
            
            MAX(ll) = LMU{ll}(N,3);
            MAX(ll) = ceil(MAX(ll)/10)*10;
            
            TYPE = zeros(N,1);
            for n = 1:N
                if isnan(DATA{ll}(m).d18O(n)) == 1
                    if isempty(DATA{ll}(m).radiocarbon{n}) == 1
                        % additional depths have label 1
                        TYPE(n) = 1;
                    else
                        % radiocarbon depths have label 2
                        TYPE(n) = 2;
                    end
                else
                    if isempty(DATA{ll}(m).radiocarbon{n}) == 1
                        % d18O depths have label 3
                        TYPE(n) = 3;
                    else
                        % depths of both kindes of data have label 4
                        TYPE(n) = 4;
                    end
                end
            end
            
            DEPTH{ll} = [Depth,TYPE];
        end
        
        DEPTH = cat(1,DEPTH{:});
        [~,order] = sort(DEPTH(:,1),'ascend');
        DEPTH = DEPTH(order,:);
        
        DEPTH_new = zeros(size(DEPTH));
        DEPTH_new(1,1) = DEPTH(1,1);
        k = 1;
        for n = 2:size(DEPTH,1)
            if abs(DEPTH(n,1)-DEPTH_new(k,1)) >= 1e-8
                k = k + 1;
                DEPTH_new(k,1) = DEPTH(n,1);
            end
        end
        DEPTH_new = DEPTH_new(1:k,:);
        for n = 1:size(DEPTH_new,1)
            index = (abs(DEPTH(:,1)-DEPTH_new(n,1))<1e-8);
            AA = DEPTH(index,2);
            if sum(AA==2) > 0 && sum(AA==3) > 0
                DEPTH_new(n,2) = 4;
            else
                DEPTH_new(n,2) = max(AA);
            end
        end
        DEPTH = DEPTH_new(:,1);
        TYPE = DEPTH_new(:,2);
        
        
        list = cell(L,1);
        h = zeros(L,1);
        fig = figure;
        
        subplot(1,3,[1,2]);
        hold on;
        tt = DATA{ll}(m).name;
        tt(tt=='_') = '-';
        title(tt,'FontSize',16);
        
        for ll = 1:L
            xx = [LMU{ll}(:,1);flipud(LMU{ll}(:,3))];
            yy = [DATA{ll}(m).depth;flipud(DATA{ll}(m).depth)];
            patch(xx,yy,1,'FaceColor',cc(ll,:),'FaceAlpha',0.1,'EdgeColor','none');
            h(ll) = plot(LMU{ll}(:,2),DATA{ll}(m).depth,'--','LineWidth',2,'Color',cc(ll,:));
            list{ll} = data_type{ll};
        end
        
        DET = 0;
        for ll = 1:L
            if isempty(CI{ll}(m).depth) == 0
                DET = 1;
                for j = 1:length(CI{ll}(m).depth)
                    h(L+1) = plot([CI{ll}(m).lower(j),CI{ll}(m).upper(j)],[CI{ll}(m).depth(j),CI{ll}(m).depth(j)],'c','LineWidth',2);
                end
            end
        end
        if DET == 1
            q = L + 1;
            list{q} = 'Individual 14C';
        else
            q = L;
        end
        
        DET = 0;
        for ll = 1:L
            for n = 1:length(SAMPLES{ll}(m).depth)
                if isnan(DATA{ll}(m).suggested_age(n,1)) == 0
                    DET = 1;
                    mu = DATA{ll}(m).suggested_age(n,1);
                    sig = DATA{ll}(m).suggested_age(n,2);
                    if DATA{ll}(m).suggested_age(n,3) == 0
                        h(q+1) = plot([max(mu-1.96*sig,0),min(mu+1.96*sig,MAX(ll))],[DATA{ll}(m).depth(n),DATA{ll}(m).depth(n)],'--m','LineWidth',2);
                    elseif DATA{ll}(m).suggested_age(n,3) == 1
                        h(q+1) = plot([max(mu-sig,0),min(mu+sig,MAX(ll))],[DATA{ll}(m).depth(n),DATA{ll}(m).depth(n)],'m','LineWidth',2);
                    end
                end
            end
        end
        if DET == 1
            q = q + 1;
            list{q} = 'Age Guess';
        end
        
        if sum(TYPE == 3) > 0
            q = q + 1;
            h(q) = plot(zeros(sum(TYPE == 3),1),DEPTH(TYPE == 3,:),'>r');
            list{q} = 'd18O';
        end
        if sum(TYPE == 2) > 0
            q = q + 1;
            h(q) = plot(zeros(sum(TYPE == 2),1),DEPTH(TYPE == 2,:),'>c');
            list{q} = '14C';
        end
        if sum(TYPE == 4) > 0
            q = q + 1;
            h(q) = plot(zeros(sum(TYPE == 4),1),DEPTH(TYPE == 4,:),'>g');
            list{q} = 'Both';
        end
        if sum(TYPE == 1) > 0
            q = q + 1;
            h(q) = plot(zeros(sum(TYPE == 1),1),DEPTH(TYPE == 1,:),'>k');
            list{q} = 'None';
        end
        
        XMAX = max(MAX);
        YMAX = -inf;
        for ll = 1:L
            YMAX = max([YMAX,DATA{ll}(m).depth(end)]);
        end
        
        xlabel('age (ky)','FontSize',12);
        ylabel('depth (m)','FontSize',12);
        xlim([0 XMAX]);
        ylim([0 YMAX*1.05]);
        legend(h,list,'Location','SouthEast','FontSize',12);
        
        %{
        MIN = -inf;
        MAX = inf;
        for ll = 1:L
            MIN = max([MIN,DATA{ll}(m).depth(1)]);
            MAX = min([MAX,DATA{ll}(m).depth(end)]);
        end
        %}
        MIN = DATA{1}(m).depth(1);
        MAX = DATA{1}(m).depth(end);
        
        DEPTH = cell(L,1);
        for ll = 1:L
            index = (DATA{ll}(m).depth>=MIN)&(DATA{ll}(m).depth<=MAX);
            DEPTH{ll} = DATA{ll}(m).depth(index);
            LMU{ll} = LMU{ll}(index,:);
        end
        
        QQ = LMU{1}(:,2);
        
        for ll = 1:L
            LMU{ll}(:,1) = LMU{ll}(:,1) - interp1(DEPTH{1},QQ,DEPTH{ll});
            LMU{ll}(:,2) = LMU{ll}(:,2) - interp1(DEPTH{1},QQ,DEPTH{ll});
            LMU{ll}(:,3) = LMU{ll}(:,3) - interp1(DEPTH{1},QQ,DEPTH{ll});
        end
        
        MAX = 6;
        MIN = -6;
        
        subplot(1,3,3);
        hold on;
        for ll = 1:L
            xx = [LMU{ll}(:,1);flipud(LMU{ll}(:,3))];
            yy = [DEPTH{ll};flipud(DEPTH{ll})];
            patch(xx,yy,1,'FaceColor',cc(ll,:),'FaceAlpha',0.1,'EdgeColor','none');
        end
        for ll = 1:L
            plot(LMU{ll}(:,2),DEPTH{ll},'--','LineWidth',2,'Color',cc(ll,:));
        end
        xlabel('relative age (ky)','FontSize',12);
        xlim([MIN MAX]);
        ylim([0 YMAX*1.05]);
        
        set(fig,'Position',[10 10 1200 500]);
        movegui(fig,'center');
    end
end


end