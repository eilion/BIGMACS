function [Samples,QQ] = init_Samples(data,core_param,global_param,cal_curve,stack,mode)

L = length(data);

Samples = struct('name',cell(L,1),'depth',cell(L,1),'ages',cell(L,1),'isoutlier',cell(L,1));
for ll = 1:L
    Samples(ll).name = data(ll).name;
    Samples(ll).depth = data(ll).depth;
end

QQ = cell(L,1);
parfor ll = 1:L
    QQ{ll} = getProposal(data(ll),core_param(ll),global_param,stack,cal_curve,mode);
end



end

