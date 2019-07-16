function AgeVsD18O(results_path)

results = load([results_path,'/results.mat']);
results = results.results;
data = results.data;
Samples = results.samples;
core_param = results.core_param;
stack = results.stack;
data_type = results.type;

L = length(data);

QQQ = cell(L,3);

for ll = 1:L
    N = length(data(ll).depth);
    
    AMAX = size(data(ll).d18O,2);
    
    H = cat(2,Samples(ll).isoutlier{:});
    H = reshape(H,[AMAX*N,size(H,2)/AMAX]);
    H = (mean(H,2)>0.5);
    H = reshape(H,[N,AMAX]);
    
    m = 0;
    OLP = cell(N,1);
    for n = 1:N
        K = sum(isnan(data(ll).d18O(n,:))==0);
        if K > 1
            m = m + 1;
            AA = quantile(Samples(ll).ages(n,:),0.5)*ones(1,K);
            BB = (data(ll).d18O(n,1:K)-core_param(ll).shift)/core_param(ll).scale;
            CC = H(n,1:K);
            [~,order] = sort(BB,'ascend');
            BB = BB(order);
            CC = CC(order);
            OLP{m} = {AA,BB,CC};
        end
    end
    OLP = OLP(1:m);
    
    M = 0;
    for n = 1:N
        if isnan(data(ll).d18O(n,1)) == 0
            M = M + 1;
        end
    end
    
    LMU = zeros(M,3);
    d18O = zeros(M,AMAX);
    IsOutlier = zeros(M,AMAX);
    m = 0;
    for n = 1:N
        if isnan(data(ll).d18O(n,1)) == 0
            m = m + 1;
            LMU(m,1) = quantile(Samples(ll).ages(n,:),0.025);
            LMU(m,2) = quantile(Samples(ll).ages(n,:),0.5);
            LMU(m,3) = quantile(Samples(ll).ages(n,:),0.975);
            d18O(m,:) = (data(ll).d18O(n,:)-core_param(ll).shift)/core_param(ll).scale;
            IsOutlier(m,:) = H(n,:);
        end
    end
    
    index = (isnan(d18O)==0);
    AAA = repmat(LMU(:,1),[1,AMAX]);
    BBB = repmat(LMU(:,2),[1,AMAX]);
    CCC = repmat(LMU(:,3),[1,AMAX]);
    
    AAA = AAA(index);
    BBB = BBB(index);
    CCC = CCC(index);
    LMU = [AAA,BBB,CCC];
    d18O = d18O(index);
    IsOutlier = IsOutlier(index);
    M = size(LMU,1);
    
    QQQ{ll,1} = LMU(:,2);
    QQQ{ll,2} = d18O;
    QQQ{ll,3} = IsOutlier;
    
    MAX = max(LMU(:,3));
    MAX = ceil(MAX/10)*10;
    
    
    if isempty(OLP) == 1
        fig = figure;
        hold on;
        tt = [data(ll).name,' (',data_type,')'];
        tt(tt=='_') = '-';
        title(tt,'FontSize',16);
        xx = [stack(:,1);flipud(stack(:,1))];
        yy1 = [stack(:,2)-stack(:,3);flipud(stack(:,2)+stack(:,3))];
        yy2 = [stack(:,2)-2*stack(:,3);flipud(stack(:,2)+2*stack(:,3))];
        patch(xx,yy1,1,'FaceColor','k','FaceAlpha',0.1,'EdgeColor','none');
        patch(xx,yy2,1,'FaceColor','k','FaceAlpha',0.1,'EdgeColor','none');
        for m = 1:M
            if IsOutlier(m) < 0.5
                plot([LMU(m,1),LMU(m,3)],[d18O(m),d18O(m)],':b','LineWidth',2);
                plot(LMU(m,2),d18O(m),'*r','LineWidth',2);
            else
                plot([LMU(m,1),LMU(m,3)],[d18O(m),d18O(m)],':r','LineWidth',2);
                plot(LMU(m,2),d18O(m),'*b','LineWidth',2);
            end
        end
        xlabel('age (ky)','FontSize',12);
        ylabel('{\delta}^{18}O (‰)','FontSize',12);
        xlim([stack(1,1) MAX]);
        ylim([2.5,5.5]);
        set(gca,'YDir','rev');
        set(fig,'Position',[10 10 1200 500]);
        movegui(fig,'center');
    else
        fig = figure;
        
        subplot(2,1,1);
        hold on;
        tt = [data(ll).name,' (',data_type,')'];
        tt(tt=='_') = '-';
        title(tt,'FontSize',16);
        xx = [stack(:,1);flipud(stack(:,1))];
        yy1 = [stack(:,2)-stack(:,3);flipud(stack(:,2)+stack(:,3))];
        yy2 = [stack(:,2)-2*stack(:,3);flipud(stack(:,2)+2*stack(:,3))];
        patch(xx,yy1,1,'FaceColor','k','FaceAlpha',0.1,'EdgeColor','none');
        patch(xx,yy2,1,'FaceColor','k','FaceAlpha',0.1,'EdgeColor','none');
        for m = 1:M
            if IsOutlier(m) < 0.5
                plot([LMU(m,1),LMU(m,3)],[d18O(m),d18O(m)],':b','LineWidth',2);
                plot(LMU(m,2),d18O(m),'*r','LineWidth',2);
            else
                plot([LMU(m,1),LMU(m,3)],[d18O(m),d18O(m)],':r','LineWidth',2);
                plot(LMU(m,2),d18O(m),'*b','LineWidth',2);
            end
        end
        xlabel('age (ky)','FontSize',12);
        ylabel('{\delta}^{18}O (‰)','FontSize',12);
        xlim([stack(1,1) MAX]);
        ylim([2.5,5.5]);
        set(gca,'YDir','rev');
        
        subplot(2,1,2);
        hold on;
        title('Multiple Observations at the Same Depths','FontSize',16);
        xx = [stack(:,1);flipud(stack(:,1))];
        yy1 = [stack(:,2)-stack(:,3);flipud(stack(:,2)+stack(:,3))];
        yy2 = [stack(:,2)-2*stack(:,3);flipud(stack(:,2)+2*stack(:,3))];
        patch(xx,yy1,1,'FaceColor','k','FaceAlpha',0.1,'EdgeColor','none');
        patch(xx,yy2,1,'FaceColor','k','FaceAlpha',0.1,'EdgeColor','none');
        for m = 1:length(OLP)
            plot(OLP{m}{1},OLP{m}{2},'--g','LineWidth',2);
            index = (OLP{m}{3}==0);
            plot(OLP{m}{1}(index),OLP{m}{2}(index),'*r','LineWidth',2);
            plot(OLP{m}{1}(~index),OLP{m}{2}(~index),'*b','LineWidth',2);
        end
        xlabel('age (ky)','FontSize',12);
        ylabel('{\delta}^{18}O (‰)','FontSize',12);
        xlim([stack(1,1) MAX]);
        ylim([2.5,5.5]);
        set(gca,'YDir','rev');
        
        set(fig,'Position',[10 10 1200 1100]);
        movegui(fig,'center');
    end
end

%{
% The below lines are not defined for general inputs.
if L > 7
    COLOR = lines(7);
    ss = (1:7)';
    tt = linspace(1,7,L)';
    cc = zeros(L,3);
    cc(:,1) = interp1(ss,COLOR(:,1),tt);
    cc(:,2) = interp1(ss,COLOR(:,2),tt);
    cc(:,3) = interp1(ss,COLOR(:,3),tt);
else
    cc = lines(L);
end

list = cell(L,1);
for ll = 1:L
    list{ll} = data(ll).name;
    list{ll}(list{ll}=='_') = '-';
end

h = zeros(L,1);
fig = figure;
hold on;
title('DNEA Stack (part)','FontSize',16);
xx = [stack(:,1);flipud(stack(:,1))];
yy1 = [stack(:,2)-stack(:,3);flipud(stack(:,2)+stack(:,3))];
yy2 = [stack(:,2)-2*stack(:,3);flipud(stack(:,2)+2*stack(:,3))];
patch(xx,yy1,1,'FaceColor','k','FaceAlpha',0.1,'EdgeColor','none');
patch(xx,yy2,1,'FaceColor','k','FaceAlpha',0.1,'EdgeColor','none');
for ll = 1:L
    h(ll) = plot(QQQ{ll,1}(QQQ{ll,3}==0),QQQ{ll,2}(QQQ{ll,3}==0),'*','Color',cc(ll,:));
    plot(QQQ{ll,1}(QQQ{ll,3}==1),QQQ{ll,2}(QQQ{ll,3}==1),'.','Color',cc(ll,:));
end
plot(stack(:,1),stack(:,2)-stack(:,3),'-.k','Linewidth',1);
plot(stack(:,1),stack(:,2)+stack(:,3),'-.k','Linewidth',1);
xlabel('age (ky)','FontSize',12);
xlim([10 50]);
ylabel('{\delta}^{18}O (‰)','FontSize',12);
ylim([3.5,5.5]);
legend(h,list,'Location','EastOutside','FontSize',12);
set(gca,'YDir','rev');
set(fig,'Position',[10 10 1200 500]);
movegui(fig,'center');
%}


end