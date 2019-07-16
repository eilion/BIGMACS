function [LOGLIK] = getLOGLIK(data,Samples,global_param,core_param,cal_curve,stack,mode)

a = 3;
b = 4;
log_tail = log(gamma(a+0.5)) - log(gamma(a)) - 0.5*log(2*pi*b);

alpha = global_param.alpha;
beta = global_param.beta;

% d = global_param.d;

L = length(data);
LOGLIK = zeros(L,1);
stack_age = stack(:,1);

for ll = 1:L
    PHI = log([core_param(ll).phi_C;core_param(ll).phi_M;core_param(ll).phi_E;core_param(ll).phi_I]);
    
    X = Samples(ll).ages;
    index = (isnan(X(:,1))==0);
    X = X(index,:);
    
    R = core_param(ll).R;
    depth = data(ll).depth(index);
    % d18O = data(ll).d18O(index) - core_param(ll).shift;
    d18O = data(ll).d18O(index,:);
    MEAN = core_param(ll).scale*stack(:,2) + core_param(ll).shift;
    STDV = core_param(ll).scale*stack(:,3);
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
    ZZ((X(2:end,:)-X(1:end-1,:))./depth_diff./RR<1.0850) = 2;
    ZZ((X(2:end,:)-X(1:end-1,:))./depth_diff./RR<0.9220) = 1;
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
        WW = (alpha-1)*log(VV) - beta*VV - log(RR(n,:));
        WW(Z(n,:)==1) = WW(Z(n,:)==1) - log(gamcdf(0.9220,alpha,1/beta));
        WW(Z(n,:)==2) = WW(Z(n,:)==2) - log(gamcdf(1.0850,alpha,1/beta)-gamcdf(0.9220,alpha,1/beta));
        WW(Z(n,:)==3) = WW(Z(n,:)==3) - log(1-gamcdf(1.0850,alpha,1/beta));
        
        LL_Table = LL_Table + WW;
    end
    
    
    if strcmp(mode,'C14') == 1
        for n = 1:N
            if isempty(C14{n}) == 0
                Table = C14{n};
                for k = 1:size(Table,1)
                    c14_mean = interp1(cal_curve{Table(k,5)}(:,1),cal_curve{Table(k,5)}(:,2),X(n,:)')';
                    c14_stdv = interp1(cal_curve{Table(k,5)}(:,1),cal_curve{Table(k,5)}(:,3),X(n,:)')';
                    LL_Table = LL_Table - (a+0.5)*log(1+(c14_mean+Table(k,3)-Table(k,1)).^2./(2*b*(c14_stdv.^2+Table(k,2)^2+Table(k,4)^2))) - 0.5*log(c14_stdv.^2+Table(k,2)^2+Table(k,4)^2) + log_tail;
                end
            end
        end
    elseif strcmp(mode,'d18O') == 1
        for m = 1:M
            index_d18O = (isnan(d18O)==0);
            
            AGE = repmat(X,[1,size(d18O,2)]);
            AGE = AGE(index_d18O);
            D18O = d18O(index_d18O);
            H = Samples(ll).isoutlier{m}(index,:);
            H = H(index_d18O);
            
            d18O_mean = interp1(stack_age,MEAN,AGE);
            d18O_stdv = interp1(stack_age,STDV,AGE);
            
            LL_Table(m) = LL_Table(m) + sum((-0.5*(D18O-d18O_mean).^2./(d18O_stdv.^2)-log(d18O_stdv)-0.5*log(2*pi)).*(1-H));
            
            %{
            AA1 = (-0.5*(D18O-d18O_mean-d*d18O_stdv).^2./(d18O_stdv.^2)-log(d18O_stdv)-0.5*log(2*pi)).*H;
            AA2 = (-0.5*(D18O-d18O_mean+d*d18O_stdv).^2./(d18O_stdv.^2)-log(d18O_stdv)-0.5*log(2*pi)).*H;
            AMAX = max([AA1,AA2],[],2);
            LL_Table(m) = LL_Table(m) + sum(AMAX-log(2)+log(exp(AA1-AMAX)+exp(AA2-AMAX)));
            %}
            
            index_age_info = (isnan(age_info(:,1))==0&(age_info(:,3)==0));
            LL_Table(m) = LL_Table(m) + sum(-0.5*(age_info(index_age_info,1)-X(index_age_info,m)).^2./(age_info(index_age_info,2).^2)-log(age_info(index_age_info,2))-0.5*log(2*pi));
        end
    elseif strcmp(mode,'both') == 1
        for n = 1:N
            if isempty(C14{n}) == 0
                Table = C14{n};
                for k = 1:size(Table,1)
                    c14_mean = interp1(cal_curve{Table(k,5)}(:,1),cal_curve{Table(k,5)}(:,2),X(n,:)')';
                    c14_stdv = interp1(cal_curve{Table(k,5)}(:,1),cal_curve{Table(k,5)}(:,3),X(n,:)')';
                    LL_Table = LL_Table - (a+0.5)*log(1+(c14_mean+Table(k,3)-Table(k,1)).^2./(2*b*(c14_stdv.^2+Table(k,2)^2+Table(k,4)^2))) - 0.5*log(c14_stdv.^2+Table(k,2)^2+Table(k,4)^2) + log_tail;
                end
            end
        end
        
        for m = 1:M
            index_d18O = (isnan(d18O)==0);
            
            AGE = repmat(X,[1,size(d18O,2)]);
            AGE = AGE(index_d18O);
            D18O = d18O(index_d18O);
            H = Samples(ll).isoutlier{m}(index,:);
            H = H(index_d18O);
            
            d18O_mean = interp1(stack_age,MEAN,AGE);
            d18O_stdv = interp1(stack_age,STDV,AGE);
            
            LL_Table(m) = LL_Table(m) + sum((-0.5*(D18O-d18O_mean).^2./(d18O_stdv.^2)-log(d18O_stdv)-0.5*log(2*pi)).*(1-H));
            
            %{
            AA1 = (-0.5*(D18O-d18O_mean-d*d18O_stdv).^2./(d18O_stdv.^2)-log(d18O_stdv)-0.5*log(2*pi)).*H;
            AA2 = (-0.5*(D18O-d18O_mean+d*d18O_stdv).^2./(d18O_stdv.^2)-log(d18O_stdv)-0.5*log(2*pi)).*H;
            AMAX = max([AA1,AA2],[],2);
            LL_Table(m) = LL_Table(m) + sum(AMAX-log(2)+log(exp(AA1-AMAX)+exp(AA2-AMAX)));
            %}
            
            index_age_info = (isnan(age_info(:,1))==0&(age_info(:,3)==0));
            LL_Table(m) = LL_Table(m) + sum(-0.5*(age_info(index_age_info,1)-X(index_age_info,m)).^2./(age_info(index_age_info,2).^2)-log(age_info(index_age_info,2))-0.5*log(2*pi));
        end
    end
    
    amax = max(LL_Table);
    LOGLIK(ll) = amax + log(sum(exp(LL_Table-amax))) - log(M);
    
    % LOGLIK(ll) = mean(LL_Table);
end


end

