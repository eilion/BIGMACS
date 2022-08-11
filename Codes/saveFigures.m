function saveFigures(inputFile,savePath)

path = ['Inputs/',inputFile,'/list_of_figures.txt'];
if exist(path,'file') ~= 2
    path = 'Defaults/list_of_figures.txt';
end
fig_list = textread(path,'%s');

if sum(strcmp(fig_list,'age_vs_d18O')) > 0
    AgeVsD18O(savePath);
end

if sum(strcmp(fig_list,'age_vs_depth')) > 0
    AgeVsDepth(savePath);
end

if sum(strcmp(fig_list,'age_vs_sedRate')) > 0
    AgeVsSedRate(savePath);
end

if sum(strcmp(fig_list,'stack_summary')) > 0
    StackSummary(savePath);
end


end