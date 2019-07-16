function [data] = getAugData(core_param,data)

L = length(core_param);

for ll = 1:L
    depth = data(ll).depth;
    d18O = data(ll).d18O;
    C14 = data(ll).radiocarbon;
    add_age = data(ll).suggested_age;
    
    N = length(depth);
    
    D = depth(end)-depth(1);
    d = D*core_param(ll).res;
    
    aug_depth = [];
    m = 0;
    for n = 1:N-1
        r = floor(abs(depth(n)-depth(n+1))/d);
        if r > 0
            A = linspace(depth(n),depth(n+1),r+2);
            aug_depth(m+1:m+r) = A(2:end-1);
            m = m + r;
        end
    end
    aug_depth = aug_depth';
    
    M = length(aug_depth);
    
    depth = [depth;aug_depth];
    d18O = [d18O;NaN*ones(M,size(d18O,2))];
    C14 = [C14;cell(M,1)];
    add_age = [add_age;NaN*ones(M,3)];
    
    [~,order] = sort(depth,'ascend');
    depth = depth(order,:);
    d18O = d18O(order,:);
    C14 = C14(order,:);
    add_age = add_age(order,:);
    
    data(ll).depth = depth;
    data(ll).d18O = d18O;
    data(ll).radiocarbon = C14;
    data(ll).suggested_age = add_age;
end


end

