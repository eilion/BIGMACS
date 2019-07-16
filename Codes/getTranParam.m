function [tran_param] = getTranParam(data,Samples,core_param)

L = length(data);
AA = zeros(4,3);

for ll = 1:L
    Age = Samples(ll).ages;
    depth = data(ll).depth;
    
    M = size(Age,2);
    
    index = (isnan(Age(:,1))==0);
    Age = Age(index,:);
    depth = depth(index);
    
    Age_diff = Age(2:end,:) - Age(1:end-1,:);
    depth_diff = depth(2:end) - depth(1:end-1);
    
    RR = zeros(size(Age(2:end,:)));
    for m = 1:size(RR,2)
        RR(:,m) = interp1(core_param(ll).R(:,1),core_param(ll).R(:,2),Age(2:end,m));
    end
    
    ZZ = Age_diff./(depth_diff.*RR);
    
    YY = 3*ones(size(ZZ));
    YY(ZZ<=1.0850) = 2;
    YY(ZZ<=0.9220) = 1;
    
    for k = 1:3
        AA(1,k) = AA(1,k) + sum(YY(end,:)==k)/M;
    end    
    for k = 1:3
        for m = 1:3
            index = (YY(2:end,:)==k&YY(1:end-1,:)==m);
            AA(k+1,m) = AA(k+1,m) + sum(sum(index))/M;
        end
    end
end

tran_param = AA./sum(AA,2);


end

