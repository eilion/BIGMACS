function [LOGLIK] = getLOGLIK(data,Samples,param,target,mode)

a = param.a_d18O;
b = param.b_d18O;
log_tail = gammaln(a+0.5) - gammaln(a) - 0.5*log(2*pi*b);

a_C14 = param.a_C14;
b_C14 = param.b_C14;

log_tail_C14 = gammaln(a_C14+0.5) - gammaln(a_C14) - 0.5*log(2*pi*b_C14);


% d = param.d;

L = length(data);
LOGLIK = zeros(L,1);
target.stack_age = target.stack(:,1);

target.stack(:,2) = target.stack(:,2)*param.scale + param.shift;
target.stack(:,3) = target.stack(:,3)*param.scale;

for ll = 1:L
    PHI = log([data(ll).phi_C;data(ll).phi_M;data(ll).phi_E;data(ll).phi_I]);
    
    X = Samples(ll).ages;
    index = (~isnan(X(:,1)));
    X = X(index,:);
    
    R = data(ll).R;
    depth = data(ll).depth(index);
    % d18O = data(ll).d18O(index) - data(ll).shift;
    d18O = data(ll).d18O(index,:);
    MEAN = data(ll).scale*target.stack(:,2) + data(ll).shift;
    STDV = data(ll).scale*target.stack(:,3);
    age_info = data(ll).suggested_age(index,:);
    depth_diff = depth(2:end) - depth(1:end-1);
    C14 = data(ll).radiocarbon(index);
    N = length(depth);
    
    M = size(X,2);
    
    LL_Table = zeros(1,M);
    
    ZZ = 3*ones(N-1,M);
    RR = zeros(N-1,M);
    for m = 1:M
        RR(:,m) = interp1(R(:,1),R(:,2),X(2:end,m));
    end
    ZZ((X(2:end,:)-X(1:end-1,:))./depth_diff./RR<1./0.9220) = 2;
    ZZ((X(2:end,:)-X(1:end-1,:))./depth_diff./RR<1./1.0850) = 1;
    ZZ((X(2:end,:)==X(1:end-1,:))) = 4;
    Z = zeros(N,M);
    Z(1:N-1,:) = ZZ;
    Z(N,:) = 4;
    
    for nn = 1:N-1
        n = N-nn;
        for m = 1:M
            LL_Table(m) = LL_Table(m) + PHI(Z(n+1,m),Z(n,m));
        end
        VV = (X(n+1,:)-X(n,:))./(RR(n,:).*depth_diff(n));
        
        % WW = interp1(data(ll).ACC_MODEL(:,1),data(ll).ACC_MODEL(:,2),VV,'linear',-56) - log(RR(n,:));
        WW = interp1(data(ll).ACC_MODEL(:,1),data(ll).ACC_MODEL(:,2),VV,'linear',-56);
        WW(Z(n,:)==1) = WW(Z(n,:)==1) - data(ll).ACC_CONTRACTION;
        WW(Z(n,:)==2) = WW(Z(n,:)==2) - data(ll).ACC_STEADY;
        WW(Z(n,:)==3) = WW(Z(n,:)==3) - data(ll).ACC_EXPANSION;
        
        LL_Table = LL_Table + WW;
    end
    
    
    if strcmp(mode,'C14') == 1
        for n = 1:N
            if isempty(C14{n}) == 0
                Table = C14{n};
                for k = 1:size(Table,1)
                    c14_mean = interp1(target.cal_curve{Table(k,5)}(:,1),target.cal_curve{Table(k,5)}(:,2),X(n,:)')';
                    c14_stdv = interp1(target.cal_curve{Table(k,5)}(:,1),target.cal_curve{Table(k,5)}(:,3),X(n,:)')';
                    
                    if Table(k,end) == 1
                        LL_Table = LL_Table - 0.5*(c14_mean+Table(k,3)-Table(k,1)).^2./(c14_stdv.^2+Table(k,2)^2+Table(k,4)^2) - 0.5*log(c14_stdv.^2+Table(k,2)^2+Table(k,4)^2) - 0.5*log(2*pi);
                    else
                        LL_Table = LL_Table - (a_C14+0.5)*log(1+(c14_mean+Table(k,3)-Table(k,1)).^2./(2*b_C14*(c14_stdv.^2+Table(k,2)^2+Table(k,4)^2))) - 0.5*log(c14_stdv.^2+Table(k,2)^2+Table(k,4)^2) + log_tail_C14;
                    end
                end
            end
        end
    elseif strcmp(mode,'d18O') == 1
        for m = 1:M
            index_d18O = (~isnan(d18O));
            
            AGE = repmat(X,[1,size(d18O,2)]);
            AGE = AGE(index_d18O);
            D18O = d18O(index_d18O);
            
            d18O_mean = interp1(target.stack_age,MEAN,AGE);
            d18O_stdv = interp1(target.stack_age,STDV,AGE);
            
            LL_Table(m) = LL_Table(m) + sum(-(a+0.5)*log(1+(D18O-d18O_mean).^2./(2*b*d18O_stdv.^2))-log(d18O_stdv)+log_tail);
            
            %{
            AA1 = (-0.5*(D18O-d18O_mean-d*d18O_stdv).^2./(d18O_stdv.^2)-log(d18O_stdv)-0.5*log(2*pi)).*H;
            AA2 = (-0.5*(D18O-d18O_mean+d*d18O_stdv).^2./(d18O_stdv.^2)-log(d18O_stdv)-0.5*log(2*pi)).*H;
            AMAX = max([AA1,AA2],[],2);
            LL_Table(m) = LL_Table(m) + sum(AMAX-log(2)+log(exp(AA1-AMAX)+exp(AA2-AMAX)));
            %}
            
            index_age_info = (~isnan(age_info(:,1))&(age_info(:,3)==0));
            LL_Table(m) = LL_Table(m) + sum(-0.5*(age_info(index_age_info,1)-X(index_age_info,m)).^2./(age_info(index_age_info,2).^2)-log(age_info(index_age_info,2))-0.5*log(2*pi));
        end
    elseif strcmp(mode,'both') == 1
        for n = 1:N
            if isempty(C14{n}) == 0
                Table = C14{n};
                for k = 1:size(Table,1)
                    c14_mean = interp1(target.cal_curve{Table(k,5)}(:,1),target.cal_curve{Table(k,5)}(:,2),X(n,:)')';
                    c14_stdv = interp1(target.cal_curve{Table(k,5)}(:,1),target.cal_curve{Table(k,5)}(:,3),X(n,:)')';
                    
                    if Table(k,end) == 1
                        LL_Table = LL_Table - 0.5*(c14_mean+Table(k,3)-Table(k,1)).^2./(c14_stdv.^2+Table(k,2)^2+Table(k,4)^2) - 0.5*log(c14_stdv.^2+Table(k,2)^2+Table(k,4)^2) - 0.5*log(2*pi);
                    else
                        LL_Table = LL_Table - (a_C14+0.5)*log(1+(c14_mean+Table(k,3)-Table(k,1)).^2./(2*b_C14*(c14_stdv.^2+Table(k,2)^2+Table(k,4)^2))) - 0.5*log(c14_stdv.^2+Table(k,2)^2+Table(k,4)^2) + log_tail_C14;
                    end
                end
            end
        end
        
        for m = 1:M
            index_d18O = (~isnan(d18O));
            
            AGE = repmat(X,[1,size(d18O,2)]);
            AGE = AGE(index_d18O);
            D18O = d18O(index_d18O);
            
            d18O_mean = interp1(target.stack_age,MEAN,AGE);
            d18O_stdv = interp1(target.stack_age,STDV,AGE);
            
            LL_Table(m) = LL_Table(m) + sum(-(a+0.5)*log(1+(D18O-d18O_mean).^2./(2*b*d18O_stdv.^2))-log(d18O_stdv)+log_tail);
            
            %{
            AA1 = (-0.5*(D18O-d18O_mean-d*d18O_stdv).^2./(d18O_stdv.^2)-log(d18O_stdv)-0.5*log(2*pi)).*H;
            AA2 = (-0.5*(D18O-d18O_mean+d*d18O_stdv).^2./(d18O_stdv.^2)-log(d18O_stdv)-0.5*log(2*pi)).*H;
            AMAX = max([AA1,AA2],[],2);
            LL_Table(m) = LL_Table(m) + sum(AMAX-log(2)+log(exp(AA1-AMAX)+exp(AA2-AMAX)));
            %}
            
            index_age_info = (~isnan(age_info(:,1))&(age_info(:,3)==0));
            LL_Table(m) = LL_Table(m) + sum(-0.5*(age_info(index_age_info,1)-X(index_age_info,m)).^2./(age_info(index_age_info,2).^2)-log(age_info(index_age_info,2))-0.5*log(2*pi));
        end
    end
    
    amax = max(LL_Table);
    LOGLIK(ll) = amax + log(sum(exp(LL_Table-amax))) - log(M);
    
    % LOGLIK(ll) = mean(LL_Table);
end


end

