function [data,param] = updateParameters(data,Samples,param,target,setting,mode,purpose)

L = length(data);

TT = param.tran_param;

if strcmp(setting.IsLearn_transition,'yes') == 1
    Grid = [1./1.0850,1./0.9220];
    phi_I = zeros(L,3);
    phi_C = zeros(L,3);
    phi_M = zeros(L,3);
    phi_E = zeros(L,3);
    for ll = 1:L
        depth = data(ll).depth;
        X = Samples(ll).ages;
        
        index = (~isnan(X(:,1)));
        depth = depth(index);
        X = X(index,:);
        
        N = length(depth);
        depth_diff = depth(2:end) - depth(1:end-1);
        
        M = size(X,2);
        X_diff = X(2:end,:) - X(1:end-1,:);
        R = zeros(size(X_diff));
        for m = 1:M
            R(:,m) = interp1(data(ll).R(:,1),data(ll).R(:,2),X(2:end,m));
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
    for ll = 1:L
        data(ll).phi_I = sum(phi_I,1) + TT(1,:);
        data(ll).phi_I = data(ll).phi_I/sum(data(ll).phi_I);
        data(ll).phi_C = sum(phi_C,1) + TT(2,:);
        data(ll).phi_C = data(ll).phi_C/sum(data(ll).phi_C);
        data(ll).phi_M = sum(phi_M,1) + TT(3,:);
        data(ll).phi_M = data(ll).phi_M/sum(data(ll).phi_M);
        data(ll).phi_E = sum(phi_E,1) + TT(4,:);
        data(ll).phi_E = data(ll).phi_E/sum(data(ll).phi_E);
    end
    for ll = 1:L
        data(ll).phi_I = data(ll).phi_I/sum(data(ll).phi_I);
        data(ll).phi_C = data(ll).phi_C/sum(data(ll).phi_C);
        data(ll).phi_M = data(ll).phi_M/sum(data(ll).phi_M);
        data(ll).phi_E = data(ll).phi_E/sum(data(ll).phi_E);
    end
end



beta1 = 0.9;
beta2 = 0.999;
epsilon = 1e-8;
gamma = 1e-3;

aa = param.a_d18O;
bb = param.b_d18O;

if strcmp(mode,'d18O') == 1 || strcmp(mode,'both') == 1
    
    L = length(data);
    AGE = target.stack(:,1);
    
    if strcmp(purpose,'stack')
        MU = target.stack(:,2);
        STDV = target.stack(:,3);
        
        Mw_G = zeros(1,2);
        Vw_G = zeros(1,2);
        
        Mw_C = zeros(L,2);
        Vw_C = zeros(L,2);
        
        SCH = zeros(10000,L);
        for ll = 1:L
            M = size(Samples(ll).ages,2);
            SCH(:,ll) = ceil(M*rand(10000,1));
        end
        
        for r = 1:10000
            PDEV_G = zeros(1,2);
            PDEV_C = zeros(L,2);
            
            GG = param.scale;
            DD = param.shift;
            
            for ll = 1:L
                m = SCH(r,ll);
                
                V = data(ll).d18O;
                A = repmat(Samples(ll).ages(:,m),[1,size(data(ll).d18O,2)]);
                ID = ~isnan(V);
                V = V(ID);
                A = A(ID);
                
                CC = data(ll).scale;
                HH = data(ll).shift;
                
                mu = interp1(AGE,MU,A);
                sig = interp1(AGE,STDV,A);
                
                QQ = (2*aa+1)*(V-CC*GG*mu-CC*DD-HH)./(2*bb*CC^2*GG^2*sig.^2+(V-CC*GG*mu-CC*DD-HH).^2);
                
                PDEV_G(1) = PDEV_G(1) + sum(QQ.*(V-CC*DD-HH)./GG-1./GG);
                PDEV_G(2) = PDEV_G(2) + sum(QQ.*CC);
                
                if strcmp(data(ll).islearn_scale,'yes') == 1
                    PDEV_C(ll,1) = sum(QQ.*(V-HH)./CC-1./CC);
                end
                if strcmp(data(ll).islearn_shift,'yes') == 1
                    PDEV_C(ll,2) = sum(QQ);
                end
            end
            
            Mw_G = beta1*Mw_G - (1-beta1).*PDEV_G;
            Vw_G = beta2*Vw_G + (1-beta2).*PDEV_G.*PDEV_G;
            
            Mw_C = beta1*Mw_C - (1-beta1).*PDEV_C;
            Vw_C = beta2*Vw_C + (1-beta2).*PDEV_C.*PDEV_C;
            
            param.scale = param.scale - (gamma*sqrt(1-beta2^r)/(1-beta1^r)).*Mw_G(1)./(sqrt(Vw_G(1))+epsilon);
            param.shift = param.shift - (gamma*sqrt(1-beta2^r)/(1-beta1^r)).*Mw_G(2)./(sqrt(Vw_G(2))+epsilon);
            
            for ll = 1:L
                data(ll).scale = data(ll).scale - (gamma*sqrt(1-beta2^r)/(1-beta1^r)).*Mw_C(ll,1)./(sqrt(Vw_C(ll,1))+epsilon);
                data(ll).shift = data(ll).shift - (gamma*sqrt(1-beta2^r)/(1-beta1^r)).*Mw_C(ll,2)./(sqrt(Vw_C(ll,2))+epsilon);
            end
        end
        
        
    elseif strcmp(purpose,'align')
        
        MU = target.stack(:,2);
        STDV = target.stack(:,3);
        
        Mw_C = zeros(L,2);
        Vw_C = zeros(L,2);
        
        SCH = zeros(10000,L);
        for ll = 1:L
            M = size(Samples(ll).ages,2);
            SCH(:,ll) = ceil(M*rand(10000,1));
        end
        
        GG = param.scale;
        DD = param.shift;
        
        for r = 1:10000
            PDEV_C = zeros(L,2);
            
            for ll = 1:L
                m = SCH(r,ll);
                
                V = data(ll).d18O;
                A = repmat(Samples(ll).ages(:,m),[1,size(data(ll).d18O,2)]);
                ID = ~isnan(V);
                V = V(ID);
                A = A(ID);
                
                CC = data(ll).scale;
                HH = data(ll).shift;
                
                mu = interp1(AGE,MU,A);
                sig = interp1(AGE,STDV,A);
                
                QQ = (2*aa+1)*(V-CC*GG*mu-CC*DD-HH)./(2*bb*CC^2*GG^2*sig.^2+(V-CC*GG*mu-CC*DD-HH).^2);
                
                if strcmp(data(ll).islearn_scale,'yes') == 1
                    PDEV_C(ll,1) = sum(QQ.*(V-HH)./CC-1./CC);
                end
                if strcmp(data(ll).islearn_shift,'yes') == 1
                    PDEV_C(ll,2) = sum(QQ);
                end
            end
            
            Mw_C = beta1*Mw_C - (1-beta1).*PDEV_C;
            Vw_C = beta2*Vw_C + (1-beta2).*PDEV_C.*PDEV_C;
            
            for ll = 1:L
                data(ll).scale = data(ll).scale - (gamma*sqrt(1-beta2^r)/(1-beta1^r)).*Mw_C(ll,1)./(sqrt(Vw_C(ll,1))+epsilon);
                data(ll).shift = data(ll).shift - (gamma*sqrt(1-beta2^r)/(1-beta1^r)).*Mw_C(ll,2)./(sqrt(Vw_C(ll,2))+epsilon);
            end
        end
    end  
end



end