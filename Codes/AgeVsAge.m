function AgeVsAge(titles,results_path)

Ages = cell(2,1);
Data = cell(2,1);
data_type = cell(2,1);

for ll = 1:2
    results = load([results_path{ll},'/results.mat']);
    results = results.results;
    Data{ll} = results.data;
    Ages{ll} = results.samples;
    data_type{ll} = results.type;
    if strcmp(data_type{ll},'both') == 1
        data_type{ll} = 'dual proxy';
    elseif strcmp(data_type{ll},'C14') == 1
        data_type{ll} = '14C';
    end
end

M = length(Ages{1});

% COLOR = jet(101);
% ss = (1:101)';
% tt = linspace(1,101,M)';
if M > 7
    COLOR = lines(7);
    ss = (1:7)';
    tt = linspace(1,7,M)';
    cc = zeros(M,3);
    cc(:,1) = interp1(ss,COLOR(:,1),tt);
    cc(:,2) = interp1(ss,COLOR(:,2),tt);
    cc(:,3) = interp1(ss,COLOR(:,3),tt);
else
    cc = lines(M);
end

MEANs = cell(M,1);
As = cell(M,1);
Bs = cell(M,1);
CXs = cell(M,1);
CYs = cell(M,1);
Angles = cell(M,1);

if strcmp(data_type{2},'14C') == 1
    MED = cell(M,1);
    for m = 1:M
        index = cellfun(@isempty,Data{2}(m).radiocarbon);
        AA = Ages{2}(m).ages(~index,:);
        MED{m} = mean(AA,2);
    end
end

DIFF = cell(M,2);

MAX = -inf;
for m = 1:M
    nSamples_1 = size(Ages{1}(m).ages,2);
    nSamples_2 = size(Ages{2}(m).ages,2);
    
    nSamples = min([nSamples_1,nSamples_2]);
    Ages{1}(m).ages = Ages{1}(m).ages(:,1:nSamples);
    Ages{2}(m).ages = Ages{2}(m).ages(:,1:nSamples);
    
    depth_1 = Data{1}(m).depth;
    depth_2 = Data{2}(m).depth;
    
    st_1 = min(depth_1);
    st_2 = min(depth_2);
    ed_1 = max(depth_1);
    ed_2 = max(depth_2);
    
    st = max(st_1,st_2);
    ed = min(ed_1,ed_2);
    
    index_1 = (depth_1>=st)&(depth_1<=ed);
    index_2 = (depth_2>=st)&(depth_2<=ed);
    
    depth_1 = depth_1(index_1);
    depth_2 = depth_2(index_2);
    age_1 = Ages{1}(m).ages(index_1,:);
    age_2 = Ages{2}(m).ages(index_2,:);
    
    AMAX1 = max(quantile(age_1(end,:),0.975));
    AMAX2 = max(quantile(age_2(end,:),0.975));
    
    MAX = max([MAX,AMAX1,AMAX2]);
    
    age_1_temp = age_2;
    for k = 1:size(age_1,2)
        age_1_temp(:,k) = interp1(depth_1,age_1(:,k),depth_2);
    end
    age_1 = age_1_temp;
    
    N = length(depth_2);
    
    %{
    age_2_temp = age_1;
    for k = 1:size(age_2,2)
        age_2_temp(:,k) = interp1(depth_2,age_2(:,k),depth_1);
    end
    age_2 = age_2_temp;
    %}
    
    % N = length(depth_1);
    
    Means = zeros(2,N);
    A = zeros(N,1);
    B = zeros(N,1);
    Cx = zeros(N,1);
    Cy = zeros(N,1);
    Angle = zeros(N,1);
    
    DIFF{m,1} = age_1 - age_2;
    DIFF{m,2} = median(age_2,2);
    
    for n = 1:N
        data_pair = [age_1(n,:);age_2(n,:)];
        Mean = sum(data_pair,2)/size(data_pair,2);
        Means(:,n) = Mean;
        COV = (data_pair-Mean)*(data_pair-Mean)'/size(data_pair,2);
        [V,D] = eig(COV);
        if D(1,1) < D(2,2)
            VV = V(:,2);
            DD = [D(2,2);D(1,1)];
        else
            VV = V(:,1);
            DD = [D(1,1);D(2,2)];
        end
        A(n) = sqrt(5.991*DD(1));
        B(n) = sqrt(5.991*DD(2));
        Cx(n) = Mean(1);
        Cy(n) = Mean(2);
        Angle(n) = atan(VV(2)/VV(1));
    end
    
    MEANs{m} = Means;
    As{m} = A;
    Bs{m} = B;
    CXs{m} = Cx;
    CYs{m} = Cy;
    Angles{m} = Angle;
end


List = cell(M,1);
for m = 1:M
    tt = Data{1}(m).name;
    tt(tt=='_') = '-';
    List{m} = tt;
end

XMAX = -inf;
XMIN = inf;
LMU = cell(M,1);
for m = 1:M
    N = size(DIFF{m,1},1);
    AA = zeros(N,3);
    for n = 1:N
        AA(n,1) = quantile(DIFF{m,1}(n,:),0.025);
        AA(n,2) = quantile(DIFF{m,1}(n,:),0.5);
        AA(n,3) = quantile(DIFF{m,1}(n,:),0.975);
    end
    LMU{m} = AA;
    XMAX = max([XMAX;max(AA(:,3))]);
    XMIN = min([XMIN;min(AA(:,1))]);
end
XMAX = max(abs(XMAX),abs(XMIN));

h = zeros(M,1);
fig = figure;
subplot(1,2,1);
hold on;
tt = [titles{1},' vs. ',titles{2}];
title(tt,'FontSize',16);
plot([0 MAX*1.05],[0 MAX*1.05],'--k','LineWidth',2);
for ll = 1:M
    Elliptic_Confidence_Region(As{ll},Bs{ll},CXs{ll},CYs{ll},Angles{ll},cc(ll,:),'lag');
end
if strcmp(data_type{2},'14C') == 1
    for ll = 1:M
        plot(zeros(length(MED{ll}),1),MED{ll},'>','Color',cc(ll,:));
    end
end
for ll = 1:M
    h(ll) = plot(MEANs{ll}(1,:),MEANs{ll}(2,:),'Linewidth',2,'Color',cc(ll,:));
end
xlim([0 MAX*1.05]);
ylim([0 MAX*1.05]);
xlabel([data_type{1},' age (ky)'],'FontSize',12);
ylabel([data_type{2},' age (ky)'],'FontSize',12);
legend(h,List,'Location','SouthEast','FontSize',12);
axis square;

subplot(1,2,2);
hold on;
title('Differences','FontSize',16);
plot([0,0],[0,MAX*1.05],'--k','LineWidth',2);
for ll = 1:M
    xx = [LMU{ll}(:,1);flipud(LMU{ll}(:,3))];
    yy = [MEANs{ll}(2,:)';flipud(MEANs{ll}(2,:)')];
    patch(xx,yy,1,'FaceColor',cc(ll,:),'FaceAlpha',0.1,'EdgeColor','none');
end
if strcmp(data_type{2},'14C') == 1
    for ll = 1:M
        plot(-XMAX*1.05*ones(length(MED{ll}),1),MED{ll},'>','Color',cc(ll,:));
    end
end
for ll = 1:M
    plot(LMU{ll}(:,2),MEANs{ll}(2,:)','LineWidth',2,'Color',cc(ll,:));
end
xlim([-XMAX*1.05 XMAX*1.05]);
ylim([0 MAX*1.05]);
xlabel([data_type{1},' age',' - ',data_type{2},' age (ky)'],'FontSize',12);
axis square;

set(fig,'Position',[20 20 1200 500]);
movegui(fig,'center');


end

