function [core_param] = Update_Param_Transition(data,Samples,core_param,global_param,stack,IsUpdate)

L = length(data);

TT = global_param.tran_param;

if strcmp(IsUpdate,'global') == 1 || strcmp(IsUpdate,'local') == 1
    Grid = [0.9220,1.0850];
    phi_I = zeros(L,3);
    phi_C = zeros(L,3);
    phi_M = zeros(L,3);
    phi_E = zeros(L,3);
    for ll = 1:L
        depth = data(ll).depth;
        X = Samples(ll).ages;
        
        index = (isnan(X(:,1))==0);
        depth = depth(index);
        X = X(index,:);
        
        N = length(depth);
        depth_diff = depth(2:end) - depth(1:end-1);
        
        M = size(X,2);
        X_diff = X(2:end,:) - X(1:end-1,:);
        R = zeros(size(X_diff));
        for m = 1:M
            R(:,m) = interp1(stack(:,1),core_param(ll).R(:,2),X(2:end,m));
        end
        X_diff_ratio = X_diff./(depth_diff.*R);
        
        Z = 2*ones(N-1,M);
        Z(X_diff_ratio < Grid(1) & X_diff_ratio > 0) = 1;
        Z(X_diff_ratio >= Grid(2)) = 3;
        
        for k = 1:3
            phi_I(ll,k) = sum(Z(end,:)==k);
        end
        for n = 1:N-2
            for k = 1:3
                phi_C(ll,k) = phi_C(ll,k) + sum(Z(n,:)==k&Z(n+1,:)==1);
                phi_M(ll,k) = phi_M(ll,k) + sum(Z(n,:)==k&Z(n+1,:)==2);
                phi_E(ll,k) = phi_E(ll,k) + sum(Z(n,:)==k&Z(n+1,:)==3);
            end
        end
        phi_I(ll,:) = phi_I(ll,:)/M;
        phi_C(ll,:) = phi_C(ll,:)/M;
        phi_M(ll,:) = phi_M(ll,:)/M;
        phi_E(ll,:) = phi_E(ll,:)/M;
    end
    if strcmp(IsUpdate,'global') == 1
        for ll = 1:L
            core_param(ll).phi_I = sum(phi_I,1) + TT(1,:);
            core_param(ll).phi_I = core_param(ll).phi_I/sum(core_param(ll).phi_I);
            core_param(ll).phi_C = sum(phi_C,1) + TT(2,:);
            core_param(ll).phi_C = core_param(ll).phi_C/sum(core_param(ll).phi_C);
            core_param(ll).phi_M = sum(phi_M,1) + TT(3,:);
            core_param(ll).phi_M = core_param(ll).phi_M/sum(core_param(ll).phi_M);
            core_param(ll).phi_E = sum(phi_E,1) + TT(4,:);
            core_param(ll).phi_E = core_param(ll).phi_E/sum(core_param(ll).phi_E);
        end
        for ll = 1:L
            core_param(ll).phi_I = core_param(ll).phi_I/sum(core_param(ll).phi_I);
            core_param(ll).phi_C = core_param(ll).phi_C/sum(core_param(ll).phi_C);
            core_param(ll).phi_M = core_param(ll).phi_M/sum(core_param(ll).phi_M);
            core_param(ll).phi_E = core_param(ll).phi_E/sum(core_param(ll).phi_E);
        end
    else
        for ll = 1:L
            core_param(ll).phi_I = phi_I(ll,:) + TT(1,:);
            core_param(ll).phi_I = core_param(ll).phi_I/sum(core_param(ll).phi_I);
            core_param(ll).phi_C = phi_C(ll,:) + TT(2,:);
            core_param(ll).phi_C = core_param(ll).phi_C/sum(core_param(ll).phi_C);
            core_param(ll).phi_M = phi_M(ll,:) + TT(3,:);
            core_param(ll).phi_M = core_param(ll).phi_M/sum(core_param(ll).phi_M);
            core_param(ll).phi_E = phi_E(ll,:) + TT(4,:);
            core_param(ll).phi_E = core_param(ll).phi_E/sum(core_param(ll).phi_E);
        end
        for ll = 1:L
            core_param(ll).phi_I = core_param(ll).phi_I/sum(core_param(ll).phi_I);
            core_param(ll).phi_C = core_param(ll).phi_C/sum(core_param(ll).phi_C);
            core_param(ll).phi_M = core_param(ll).phi_M/sum(core_param(ll).phi_M);
            core_param(ll).phi_E = core_param(ll).phi_E/sum(core_param(ll).phi_E);
        end
    end
end


end

