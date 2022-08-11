function StackSummary(results_path)

results = load([results_path,'/results.mat']);

data = results.data;
Samples = results.samples;
stack = results.target.stack;
old_stack = results.target.init_stack;
setting = results.setting;

data_type = setting.data_type;

if ~strcmp(data_type,'C14')
    
    
    path = [results_path,'/figures/stack_summary'];
    mkdir(path);
    
    L = length(data);
    
    list = cell(L,1);
    for ll = 1:L
        list{ll} = data(ll).name;
        list{ll}(list{ll}=='_') = '-';
    end
    
    COLOR = jet(101);
    ss = (1:101)';
    qq = linspace(1,101,L)';
    cc = zeros(L,3);
    cc(:,1) = interp1(ss,COLOR(:,1),qq);
    cc(:,2) = interp1(ss,COLOR(:,2),qq);
    cc(:,3) = interp1(ss,COLOR(:,3),qq);
    
    
    % Plot Stack vs. Initial Stack:
    h = zeros(2+L,1);
    fig = figure;
    hold on;
    title('New Stack vs. Initial Stack','FontSize',16);
    xx = [old_stack(:,1);flipud(old_stack(:,1))];
    yy = [old_stack(:,2)-1.96*old_stack(:,3);flipud(old_stack(:,2)+1.96*old_stack(:,3))];
    patch(xx,yy,1,'FaceColor','g','FaceAlpha',0.1,'EdgeColor','none');
    
    for ll = 1:L
        for k = 1:size(data(ll).d18O,2)
            h(2+ll) = plot(Samples(ll).ages(:,k),data(ll).d18O(:,k),'*','Color',cc(ll,:));
        end
    end
    
    h(2) = plot(old_stack(:,1),old_stack(:,2),'--g','LineWidth',2);
    h(1) = plot(stack(:,1),stack(:,2),'--k','LineWidth',2);
    plot(stack(:,1),stack(:,2)-1.96*stack(:,3),':k','LineWidth',2);
    plot(stack(:,1),stack(:,2)+1.96*stack(:,3),':k','LineWidth',2);
    
    xlabel('age (ky)','FontSize',12);
    ylabel('{\delta}^{18}O (‰)','FontSize',12);
    xlim([0 stack(end,1)]);
    legend(h,[{'new stack';'initial stack'};list],'Location','EastOutside');
    set(fig,'Position',[10 10 1200 500]);
    set(gca,'YDir','rev');
    movegui(fig,'center');
    
    path = [results_path,'/figures/stack_summary/new_stack_vs_initial_stack.fig'];
    savefig(fig,path);
    
    
    
    % Plot the stack with all data:
    AGE = cell(L,1);
    D18O = cell(L,1);
    for ll = 1:L
        age = median(Samples(ll).ages,2);
        AMAX = size(data(ll).d18O,2);
        age = repmat(age,[1,AMAX]);
        index = (~isnan(data(ll).d18O));
        AGE{ll} = age(index);
        D18O{ll} = (data(ll).d18O(index)-data(ll).shift)/data(ll).scale;
    end
    h = zeros(L,1);
    fig = figure;
    hold on;
    xx = [stack(:,1);flipud(stack(:,1))];
    yy = [stack(:,2)-1.96*stack(:,3);flipud(stack(:,2)+1.96*stack(:,3))];
    patch(xx,yy,1,'FaceColor','g','FaceAlpha',0.15,'EdgeColor','none');
    plot(stack(:,1),stack(:,2)-1.96*stack(:,3),':k','LineWidth',1);
    plot(stack(:,1),stack(:,2)+1.96*stack(:,3),':k','LineWidth',1);
    for ll = 1:L
        h(ll) = plot(AGE{ll},D18O{ll},'*','Color',cc(ll,:),'LineWidth',2);
    end
    xlabel('age (ky)','FontSize',12);
    ylabel('{\delta}^{18}O (‰)','FontSize',12);
    xlim([0 stack(end,1)]);
    legend(h,list,'Location','NorthEast');
    set(fig,'Position',[10 10 800 600]);
    set(gca,'YDir','rev');
    movegui(fig,'center');
    
    path = [results_path,'/figures/stack_summary/stack_and_d18O.fig'];
    savefig(fig,path);
    
    
    
    % Plot the histogram of normalized d18Os:
    AGE = cell(L,1);
    D18O = cell(L,1);
    M = inf;
    for ll = 1:L
        M = min([M,size(Samples(ll).ages,2)]);
    end
    
    for ll = 1:L
        AMAX = size(data(ll).d18O,2);
        
        Age = Samples(ll).ages(:,1:M);
        d18O = (data(ll).d18O-data(ll).shift)/data(ll).scale;
        index = (~isnan(d18O));
        
        XX = cell(M,1);
        YY = cell(M,1);
        for m = 1:M
            XX{m} = repmat(Age(:,m),[1,AMAX]);
            XX{m} = XX{m}(index);
            YY{m} = d18O(index);
        end
        AGE{ll} = cat(1,XX{:});
        D18O{ll} = cat(1,YY{:});
    end
    X = cat(1,AGE{:});
    Y = cat(1,D18O{:});
    mu = interp1(stack(:,1),stack(:,2),X);
    sig = interp1(stack(:,1),stack(:,3),X);
    Y = (Y-mu)./sig;
    
    Z = cat(1,D18O{:});
    mu = interp1(old_stack(:,1),old_stack(:,2),X);
    sig = interp1(old_stack(:,1),old_stack(:,3),X);
    Z = (Z-mu)./sig;
    
    h = zeros(3,1);
    fig = figure;
    hold on;
    title('Histograms of Standardized {\delta}^{18}O Residuals','FontSize',16);
    bins = -6:0.25:6;
    x = -6:0.01:6;
    y = normpdf(x,0,1);
    hh2 = histogram(Z,bins,'FaceColor','c');
    hh1 = histogram(Y,bins,'FaceColor','m');
    hh1.Normalization = 'pdf';
    hh2.Normalization = 'pdf';
    h(2) = hh1;
    h(1) = plot(x,y,'k','LineWidth',3);
    h(3) = hh2;
    xlim([-6 6]);
    ylim([0 0.75]);
    legend(h,{'Std. Normal','new stack','initial stack'},'Location','EastOutside','FontSize',12);
    legend(h,{'Std. Normal','DNEA','DNA'},'Location','EastOutside','FontSize',12);
    set(fig,'Position',[10 10 1200 500]);
    movegui(fig,'center');
    
    path = [results_path,'/figures/stack_summary/residual_histogram.fig'];
    savefig(fig,path);
    
end

close all;


end