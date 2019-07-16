function Stack_Comparison(Titles,results_path)

L = length(results_path);

COLOR = jet(101);
ss = (1:101)';
qq = linspace(1,101,L)';
cc = zeros(L,3);
cc(:,1) = interp1(ss,COLOR(:,1),qq);
cc(:,2) = interp1(ss,COLOR(:,2),qq);
cc(:,3) = interp1(ss,COLOR(:,3),qq);


MIN = inf;
MAX = -inf;
Stacks = cell(L,1);
for ll = 1:L
    path = [results_path{ll},'/stack.txt'];
    Stacks{ll} = getStack(path);
    MIN = min(MIN,min(Stacks{ll}(:,1)));
    MAX = max(MAX,max(Stacks{ll}(:,1)));
    Titles{ll}(Titles{ll}=='_') = '-';
end

h = zeros(L,1);
fig = figure;
hold on;
title('Stack Comparison');
for ll = 1:L
    plot(Stacks{ll}(:,1),Stacks{ll}(:,2)-1.96*Stacks{ll}(:,3),':','Color',cc(ll,:),'LineWidth',2);
    plot(Stacks{ll}(:,1),Stacks{ll}(:,2)+1.96*Stacks{ll}(:,3),':','Color',cc(ll,:),'LineWidth',2);
    h(ll) = plot(Stacks{ll}(:,1),Stacks{ll}(:,2),'--','Color',cc(ll,:),'LineWidth',2);
end
xlabel('age (ky)','FontSize',12);
ylabel('{\delta}^{18}O (?)','FontSize',12);
xlim([MIN MAX]);
ylim([2.5 5.5]);
legend(h,Titles,'Location','EastOutside');
set(fig,'Position',[10 10 1200 500]);
set(gca,'YDir','rev');
movegui(fig,'center');


end