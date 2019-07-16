function Stack_Summary(results_path)

results = load([results_path,'/results.mat']);
results = results.results;
data = results.data;
Samples = results.samples;
core_param = results.core_param;
stack = results.stack;
old_stack = results.old_stack;
rcd_stack = results.rcd_stack;

tt = results.name;
tt(tt=='_') = '-';

L = length(rcd_stack);

COLOR = jet(101);
ss = (1:101)';
qq = linspace(1,101,L)';
cc = zeros(L,3);
cc(:,1) = interp1(ss,COLOR(:,1),qq);
cc(:,2) = interp1(ss,COLOR(:,2),qq);
cc(:,3) = interp1(ss,COLOR(:,3),qq);


% Plot Stack vs. Initial Stack:
h = zeros(2,1);
fig = figure;
hold on;
title([tt,' vs. Initial Stack'],'FontSize',16);
xx = [old_stack(:,1);flipud(old_stack(:,1))];
yy = [old_stack(:,2)-1.96*old_stack(:,3);flipud(old_stack(:,2)+1.96*old_stack(:,3))];
patch(xx,yy,1,'FaceColor','g','FaceAlpha',0.1,'EdgeColor','none');
% xx = [stack(:,1);flipud(stack(:,1))];
% yy = [stack(:,2)-1.96*stack(:,3);flipud(stack(:,2)+1.96*stack(:,3))];
% patch(xx,yy,1,'FaceColor','k','FaceAlpha',0.1,'EdgeColor','none');
h(2) = plot(old_stack(:,1),old_stack(:,2),'--g','LineWidth',2);
h(1) = plot(stack(:,1),stack(:,2),'--k','LineWidth',2);
plot(stack(:,1),stack(:,2)-1.96*stack(:,3),':k','LineWidth',2);
plot(stack(:,1),stack(:,2)+1.96*stack(:,3),':k','LineWidth',2);
xlabel('age (ky)','FontSize',12);
ylabel('{\delta}^{18}O (‰)','FontSize',12);
xlim([0 stack(end,1)]);
ylim([2.5 5.5]);
legend(h,{tt,'initial stack'},'Location','EastOutside');
set(fig,'Position',[10 10 1200 500]);
set(gca,'YDir','rev');
movegui(fig,'center');


% Plot Stack vs. Record-specific Stacks:
list = cell(L+1,1);
list{1} = 'local stack';
for ll = 1:L
    list{ll+1} = rcd_stack(ll).name;
    list{ll+1}(list{ll+1}=='_') = '-';
end
h = zeros(L+1,1);
fig = figure;
hold on;
title([tt,' vs. Record-specific GP Models'],'FontSize',16);
for ll = 1:L
    xx = [rcd_stack(ll).age;flipud(rcd_stack(ll).age)];
    yy = [rcd_stack(ll).mean-1.96*rcd_stack(ll).stdv;flipud(rcd_stack(ll).mean+1.96*rcd_stack(ll).stdv)];
    patch(xx,yy,1,'FaceColor',cc(ll,:),'FaceAlpha',0.1,'EdgeColor','none');
end
% xx = [stack(:,1);flipud(stack(:,1))];
% yy = [stack(:,2)-1.96*stack(:,3);flipud(stack(:,2)+1.96*stack(:,3))];
% patch(xx,yy,1,'FaceColor','k','FaceAlpha',0.1,'EdgeColor','none');
for ll = 1:L
    h(ll+1) = plot(rcd_stack(ll).age,rcd_stack(ll).mean,'--','LineWidth',1,'Color',cc(ll,:));
end
h(1) = plot(stack(:,1),stack(:,2),'--k','LineWidth',2);
plot(stack(:,1),stack(:,2)-1.96*stack(:,3),':k','LineWidth',2);
plot(stack(:,1),stack(:,2)+1.96*stack(:,3),':k','LineWidth',2);
xlabel('age (ky)','FontSize',12);
ylabel('{\delta}^{18}O (‰)','FontSize',12);
xlim([0 stack(end,1)]);
ylim([2.5 5.5]);
legend(h,list,'Location','EastOutside');
set(fig,'Position',[10 10 1200 500]);
set(gca,'YDir','rev');
movegui(fig,'center');

for ll = 1:L
    N = length(data(ll).depth);
    
    M = 0;
    for n = 1:N
        if isnan(data(ll).d18O(n,1)) == 0
            M = M + 1;
        end
    end
    AMAX = size(data(ll).d18O,2);
    
    H = cat(2,Samples(ll).isoutlier{:});
    H = reshape(H,[AMAX*N,size(H,2)/AMAX]);
    H = (mean(H,2)>0.5);
    H = reshape(H,[N,AMAX]);
    
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
    
    
    fig = figure;
    
    subplot(2,1,1);
    hold on;
    ss = rcd_stack(ll).name;
    ss(ss=='_') = '-';
    title([tt,' vs. ',ss,' GP Model'],'FontSize',16);
    xx = [rcd_stack(ll).age;flipud(rcd_stack(ll).age)];
    yy = [rcd_stack(ll).mean-1.96*rcd_stack(ll).stdv;flipud(rcd_stack(ll).mean+1.96*rcd_stack(ll).stdv)];
    patch(xx,yy,1,'FaceColor',cc(ll,:),'FaceAlpha',0.1,'EdgeColor','none');
    plot(rcd_stack(ll).age,rcd_stack(ll).mean,'--','LineWidth',2,'Color',cc(ll,:));
    plot(stack(:,1),stack(:,2),'--k','LineWidth',2);
    plot(stack(:,1),stack(:,2)-1.96*stack(:,3),':k','LineWidth',2);
    plot(stack(:,1),stack(:,2)+1.96*stack(:,3),':k','LineWidth',2);
    xlabel('age (ky)','FontSize',12);
    ylabel('{\delta}^{18}O (‰)','FontSize',12);
    AMIN = max([stack(1,1),floor(rcd_stack(ll).age(1)/10)*10]);
    AMAX = min([stack(end,1),ceil(rcd_stack(ll).age(end)/10)*10]);
    xlim([AMIN AMAX]);
    ylim([2.5 5.5]);
    set(gca,'YDir','rev');
    
    subplot(2,1,2);
    hold on;
    vv = [data(ll).name,' on its Record-specific GP Model'];
    vv(vv=='_') = '-';
    title(vv,'FontSize',16);
    xx = [rcd_stack(ll).age;flipud(rcd_stack(ll).age)];
    yy = [rcd_stack(ll).mean-1.96*rcd_stack(ll).stdv;flipud(rcd_stack(ll).mean+1.96*rcd_stack(ll).stdv)];
    patch(xx,yy,1,'FaceColor',cc(ll,:),'FaceAlpha',0.1,'EdgeColor','none');
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
    AMIN = max([stack(1,1),floor(rcd_stack(ll).age(1)/10)*10]);
    AMAX = min([stack(end,1),ceil(rcd_stack(ll).age(end)/10)*10]);
    xlim([AMIN AMAX]);
    ylim([2.5 5.5]);
    
    
    set(fig,'Position',[10 10 1200 1100]);
    set(gca,'YDir','rev');
    movegui(fig,'center');
end


% Plot the stack with all data:
AGE = cell(L,1);
D18O = cell(L,1);
for ll = 1:L
    age = median(Samples(ll).ages,2);
    AMAX = size(data(ll).d18O,2);
    age = repmat(age,[1,AMAX]);
    index = (isnan(data(ll).d18O)==0);
    AGE{ll} = age(index);
    D18O{ll} = (data(ll).d18O(index)-core_param(ll).shift)/core_param(ll).scale;
end
h = zeros(L+1,1);
fig = figure;
hold on;
title([tt,' with d18O Data'],'FontSize',16);
xx = [stack(:,1);flipud(stack(:,1))];
yy = [stack(:,2)-1.96*stack(:,3);flipud(stack(:,2)+1.96*stack(:,3))];
patch(xx,yy,1,'FaceColor','k','FaceAlpha',0.1,'EdgeColor','none');
h(1) = plot(stack(:,1),stack(:,2)-1.96*stack(:,3),':k','LineWidth',1);
plot(stack(:,1),stack(:,2)+1.96*stack(:,3),':k','LineWidth',1);
for ll = 1:L
    h(ll+1) = plot(AGE{ll},D18O{ll},'.','Color',cc(ll,:));
end
xlabel('age (ky)','FontSize',12);
ylabel('{\delta}^{18}O (‰)','FontSize',12);
xlim([0 stack(end,1)]);
ylim([2.5 5.5]);
legend(h,list,'Location','EastOutside');
set(fig,'Position',[10 10 1200 500]);
set(gca,'YDir','rev');
movegui(fig,'center');


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
    d18O = (data(ll).d18O-core_param(ll).shift)/core_param(ll).scale;
    index = (isnan(d18O)==0);
    
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
title('Histograms of Normalized d18O Data','FontSize',16);
bins = -6:0.25:6;
x = -6:0.01:6;
y = normpdf(x,0,1);
h(1) = plot(x,y,'g','LineWidth',3);
hh2 = histogram(Z,bins,'FaceColor','c');
hh1 = histogram(Y,bins,'FaceColor','m');
hh1.Normalization = 'pdf';
hh2.Normalization = 'pdf';
h(2) = hh1;
h(3) = hh2;
xlim([-6 6]);
ylim([0 1]);
legend(h,{'Standard Normal','new stack','initial stack'},'Location','NorthEast','FontSize',12);
set(fig,'Position',[10 10 1200 500]);
movegui(fig,'center');


end