function [Samples,QQ] = MCMC_MH(data,Samples,param,target,mode,nsamples)

max_iters = param.max_iters;

HH = param.shift;
CC = param.scale;

target.stack(:,2) = target.stack(:,2)*CC + HH;
target.stack(:,3) = target.stack(:,3)*CC;

a = param.a_d18O;
b = param.b_d18O;

a_C14 = param.a_C14;
b_C14 = param.b_C14;

% log_tail = log(gamma(a+0.5)) - log(gamma(a)) - 0.5*log(2*pi*b);

S = param.nParticles;

q = param.q;
d = param.d;

phi_I = data.phi_I;
phi_C = data.phi_C;
phi_M = data.phi_M;
phi_E = data.phi_E;
PHI = [phi_C;phi_M;phi_E;phi_I];
PHI = log(PHI);

depth = data.depth;
depth_diff = depth(2:end) - depth(1:end-1);
% d18O = (data.d18O-data.shift)/data.scale;
d18O = data.d18O;
target.stack(:,2) = data.scale*target.stack(:,2) + data.shift;
target.stack(:,3) = data.scale*target.stack(:,3);
C14 = data.radiocarbon;
Age_Info = data.suggested_age;

MAX = size(d18O,2);

R = data.R;

N = length(depth);

aa = 1:N;
index_0_next = (rem(aa,3)==1&aa>3);
index_0 = (rem(aa,3)==0&aa>2&aa<N);
index_0_prev = (rem(aa,3)==2&aa>1&aa<N-1);
index_0_prev_prev = (rem(aa,3)==1&aa<N-2);
index_1_next = (rem(aa,3)==2&aa>3);
index_1 = (rem(aa,3)==1&aa>2&aa<N);
index_1_prev = (rem(aa,3)==0&aa>1&aa<N-1);
index_1_prev_prev = (rem(aa,3)==2&aa<N-2);
index_2_next = (rem(aa,3)==0&aa>3);
index_2 = (rem(aa,3)==2&aa>2&aa<N);
index_2_prev = (rem(aa,3)==1&aa>1&aa<N-1);
index_2_prev_prev = (rem(aa,3)==0&aa<N-2);
N0 = sum(index_0);
N1 = sum(index_1);
N2 = sum(index_2);

Age_Info_0 = Age_Info(index_0,:);
Age_Info_1 = Age_Info(index_1,:);
Age_Info_2 = Age_Info(index_2,:);

A = Samples.ages;
M = size(A,2);

if nsamples > M
    index = ceil(M*rand(1,nsamples));
    A = A(:,index);
end
M = size(A,2);

ZZ = 3*ones(N-1,M);
RR = zeros(N-1,M);
for m = 1:M
    RR(:,m) = interp1(R(:,1),R(:,2),A(2:end,m));
end
ZZ((A(2:end,:)-A(1:end-1,:))./(depth_diff.*RR)<=1./0.9220) = 2;
ZZ((A(2:end,:)-A(1:end-1,:))./(depth_diff.*RR)<=1./1.0850) = 1;
Z = zeros(N,M);
Z(1:N-1,:) = ZZ;
Z(N,:) = 4;

for iters = 1:max_iters
    Rand_Seed = log(rand(N,M));
    
    % n == 1:
    A_old = A(1,:);
    A_new = normrnd(A(1,:),(A(2,:)-A(1,:)+1)/8);
    
    Z1_old = Z(1,:);
    Z1_new = 3*ones(1,M);
    RR = interp1(R(:,1),R(:,2),A(2,:));
    Z1_new((A(2,:)-A_new)./(depth_diff(1).*RR)<=1./0.9220) = 2;
    Z1_new((A(2,:)-A_new)./(depth_diff(1).*RR)<=1./1.0850) = 1;
    
    LL_old = zeros(1,M);
    LL_new = zeros(1,M);
    
    for k = 1:3
        for m = 1:3
            LL_old(Z(2,:)==k&Z1_old==m) = PHI(k,m);
            LL_new(Z(2,:)==k&Z1_new==m) = PHI(k,m);
        end
    end
    
    LL_new(A_new>=A(2,:)|A_new<target.stack(1,1)) = -inf;
    
    TT_old = (A(2,:)-A_old)./(RR.*depth_diff(1));
    TT_new = (A(2,:)-A_new)./(RR.*depth_diff(1));
    
    LL_old(TT_old>0) = LL_old(TT_old>0) + interp1(data.ACC_MODEL(:,1),data.ACC_MODEL(:,2),TT_old(TT_old>0),'linear',-56);
    LL_new(TT_new>0) = LL_new(TT_new>0) + interp1(data.ACC_MODEL(:,1),data.ACC_MODEL(:,2),TT_new(TT_new>0),'linear',-56);
    
    LL_old(TT_old<=0) = -inf;
    LL_new(TT_new<=0) = -inf;
    
    LL_old(TT_old>1./data.lower_sedrate|TT_old<1./data.upper_sedrate) = -inf;
    LL_new(TT_new>1./data.lower_sedrate|TT_new<1./data.upper_sedrate) = -inf;
    
    LL_new = LL_new - 0.5*(A_new-A_old).^2./((A(2,:)-A_new+1).^2/64) - log((A(2,:)-A_new+1)/8);
    LL_old = LL_old - 0.5*(A_new-A_old).^2./((A(2,:)-A_old+1).^2/64) - log((A(2,:)-A_old+1)/8);
    
    LL_old(Z1_old==1) = LL_old(Z1_old==1) - data.ACC_CONTRACTION;
    LL_old(Z1_old==2) = LL_old(Z1_old==2) - data.ACC_STEADY;
    LL_old(Z1_old==3) = LL_old(Z1_old==3) - data.ACC_EXPANSION;
    
    LL_new(Z1_new==1) = LL_new(Z1_new==1) - data.ACC_CONTRACTION;
    LL_new(Z1_new==2) = LL_new(Z1_new==2) - data.ACC_STEADY;
    LL_new(Z1_new==3) = LL_new(Z1_new==3) - data.ACC_EXPANSION;
    
    if ~strcmp(mode,'C14') && ~isnan(d18O(1,1))
        d18O_mu = interp1(target.stack(:,1),target.stack(:,2),A_old);
        d18O_stdv = interp1(target.stack(:,1),target.stack(:,3),A_old);
        NN = sum(~isnan(d18O(1,:)));
        for n = 1:NN
            % LL_d18O = (1-q)*exp(-(d18O(1,n)-d18O_mu).^2./(2*d18O_stdv.^2))./sqrt(2*pi*d18O_stdv.^2) + q/2*exp(-(d18O(1,n)-d18O_mu-d*d18O_stdv).^2./(2*d18O_stdv.^2))./sqrt(2*pi*d18O_stdv.^2) + q/2*exp(-(d18O(1,n)-d18O_mu+d*d18O_stdv).^2./(2*d18O_stdv.^2))./sqrt(2*pi*d18O_stdv.^2);
            % LL_old = LL_old + log(LL_d18O);
            LL_old = LL_old - (a+0.5)*log(1+(d18O(1,n)-d18O_mu).^2./(2*b*d18O_stdv.^2)) - log(d18O_stdv);
        end
        
        d18O_mu = interp1(target.stack(:,1),target.stack(:,2),A_new);
        d18O_stdv = interp1(target.stack(:,1),target.stack(:,3),A_new);
        for n = 1:NN
            % LL_d18O = (1-q)*exp(-(d18O(1,n)-d18O_mu).^2./(2*d18O_stdv.^2))./sqrt(2*pi*d18O_stdv.^2) + q/2*exp(-(d18O(1,n)-d18O_mu-d*d18O_stdv).^2./(2*d18O_stdv.^2))./sqrt(2*pi*d18O_stdv.^2) + q/2*exp(-(d18O(1,n)-d18O_mu+d*d18O_stdv).^2./(2*d18O_stdv.^2))./sqrt(2*pi*d18O_stdv.^2);
            % LL_new = LL_new + log(LL_d18O);
            LL_new = LL_new - (a+0.5)*log(1+(d18O(1,n)-d18O_mu).^2./(2*b*d18O_stdv.^2)) - log(d18O_stdv);
        end
    end
    
    if ~strcmp(mode,'d18O') && ~isempty(C14{1})
        for ll = 1:size(C14{1},1)
            c14_mu = interp1(target.cal_curve{C14{1}(ll,5)}(:,1),target.cal_curve{C14{1}(ll,5)}(:,2),A_old);
            c14_stdv = interp1(target.cal_curve{C14{1}(ll,5)}(:,1),target.cal_curve{C14{1}(ll,5)}(:,3),A_old);
            if C14{1}(ll,end) == 1
                LL_old = LL_old - 0.5*(c14_mu+C14{1}(ll,3)-C14{1}(ll,1)).^2./(c14_stdv.^2+C14{1}(ll,2)^2+C14{1}(ll,4)^2) - 0.5*log(c14_stdv.^2+C14{1}(ll,2)^2+C14{1}(ll,4)^2) - 0.5*log(2*pi);
            else
                LL_old = LL_old - (a_C14+0.5)*log(1+(c14_mu+C14{1}(ll,3)-C14{1}(ll,1)).^2./(2*b_C14*(c14_stdv.^2+C14{1}(ll,2)^2+C14{1}(ll,4)^2))) - 0.5*log(c14_stdv.^2+C14{1}(ll,2)^2+C14{1}(ll,4)^2);
            end
            
            c14_mu = interp1(target.cal_curve{C14{1}(ll,5)}(:,1),target.cal_curve{C14{1}(ll,5)}(:,2),A_new);
            c14_stdv = interp1(target.cal_curve{C14{1}(ll,5)}(:,1),target.cal_curve{C14{1}(ll,5)}(:,3),A_new);
            if C14{1}(ll,end) == 1
                LL_new = LL_new - 0.5*(c14_mu+C14{1}(ll,3)-C14{1}(ll,1)).^2./(c14_stdv.^2+C14{1}(ll,2)^2+C14{1}(ll,4)^2) - 0.5*log(c14_stdv.^2+C14{1}(ll,2)^2+C14{1}(ll,4)^2) - 0.5*log(2*pi);
            else
                LL_new = LL_new - (a_C14+0.5)*log(1+(c14_mu+C14{1}(ll,3)-C14{1}(ll,1)).^2./(2*b_C14*(c14_stdv.^2+C14{1}(ll,2)^2+C14{1}(ll,4)^2))) - 0.5*log(c14_stdv.^2+C14{1}(ll,2)^2+C14{1}(ll,4)^2);
            end
        end
    end
    
    if ~isnan(Age_Info(1,1))
        if Age_Info(1,3) == 0
            LL_old = LL_old - (A_old-Age_Info(1,1)).^2./(2*Age_Info(1,2).^2);
            LL_new = LL_new - (A_new-Age_Info(1,1)).^2./(2*Age_Info(1,2).^2);
        elseif Age_Info(1,3) == 1
            index_old = (A_old>Age_Info(1,1)-Age_Info(1,2))&(A_old<Age_Info(1,1)+Age_Info(1,2));
            index_new = (A_new>Age_Info(1,1)-Age_Info(1,2))&(A_new<Age_Info(1,1)+Age_Info(1,2));
            LL_old(~index_old) = -inf;
            LL_new(~index_new) = -inf;
        end
    end
    
    if ~isnan(data.max)
        LL_old(A_old>data.max) = -inf;
        LL_new(A_new>data.max) = -inf;
    end
    if ~isnan(data.min)
        LL_old(A_old<data.min) = -inf;
        LL_new(A_new<data.min) = -inf;
    end
    
    index = (LL_new-LL_old>Rand_Seed(1,:));
    A(1,index) = A_new(index);
    Z(1,index) = Z1_new(index);
    
    
    % n == 2:
    A_old = A(2,:);
    A_new = normrnd(A(2,:),(A(3,:)-A(1,:))/8);
    
    Z1_old = Z(1,:);
    Z1_new = 3*ones(1,M);
    RR = interp1(R(:,1),R(:,2),A_new);
    Z1_new((A_new-A(1,:))./(depth_diff(1).*RR)<=1./0.9220) = 2;
    Z1_new((A_new-A(1,:))./(depth_diff(1).*RR)<=1./1.0850) = 1;
    
    Z2_old = Z(2,:);
    Z2_new = 3*ones(1,M);
    RR = interp1(R(:,1),R(:,2),A(3,:));
    Z2_new((A(3,:)-A_new)./(depth_diff(2).*RR)<=1./0.9220) = 2;
    Z2_new((A(3,:)-A_new)./(depth_diff(2).*RR)<=1./1.0850) = 1;
    
    LL_old = zeros(1,M);
    LL_new = zeros(1,M);
    
    for k = 1:3
        for m = 1:3
            LL_old(Z(3,:)==k&Z2_old==m) = PHI(k,m);
            LL_new(Z(3,:)==k&Z2_new==m) = PHI(k,m);
        end
    end
    
    for k = 1:3
        for m = 1:3
            LL_old(Z2_old==k&Z1_old==m) = LL_old(Z2_old==k&Z1_old==m) + PHI(k,m);
            LL_new(Z2_new==k&Z1_new==m) = LL_new(Z2_new==k&Z1_new==m) + PHI(k,m);
        end
    end
    
    LL_new(A_new>=A(3,:)|A_new<=A(1,:)) = -inf;
    
    TT_old = (A(3,:)-A_old)./(RR.*depth_diff(2));
    TT_new = (A(3,:)-A_new)./(RR.*depth_diff(2));
    
    LL_old(TT_old>0) = LL_old(TT_old>0) + interp1(data.ACC_MODEL(:,1),data.ACC_MODEL(:,2),TT_old(TT_old>0),'linear',-56);
    LL_new(TT_new>0) = LL_new(TT_new>0) + interp1(data.ACC_MODEL(:,1),data.ACC_MODEL(:,2),TT_new(TT_new>0),'linear',-56);
    
    LL_old(TT_old<=0) = -inf;
    LL_new(TT_new<=0) = -inf;
    
    LL_old(TT_old>1./data.lower_sedrate|TT_old<1./data.upper_sedrate) = -inf;
    LL_new(TT_new>1./data.lower_sedrate|TT_new<1./data.upper_sedrate) = -inf;
    
    RR_old = interp1(R(:,1),R(:,2),A_old);
    RR_new = interp1(R(:,1),R(:,2),A_new);
    
    TT_old = (A_old-A(1,:))./(RR_old.*depth_diff(1));
    TT_new = (A_new-A(1,:))./(RR_new.*depth_diff(1));
    
    LL_old(TT_old>0) = LL_old(TT_old>0) + interp1(data.ACC_MODEL(:,1),data.ACC_MODEL(:,2),TT_old(TT_old>0),'linear',-56);
    LL_new(TT_new>0) = LL_new(TT_new>0) + interp1(data.ACC_MODEL(:,1),data.ACC_MODEL(:,2),TT_new(TT_new>0),'linear',-56);
    
    LL_old(TT_old<=0) = -inf;
    LL_new(TT_new<=0) = -inf;
    
    LL_old(TT_old>1./data.lower_sedrate|TT_old<1./data.upper_sedrate) = -inf;
    LL_new(TT_new>1./data.lower_sedrate|TT_new<1./data.upper_sedrate) = -inf;
    
    LL_old(Z2_old==1) = LL_old(Z2_old==1) - data.ACC_CONTRACTION;
    LL_old(Z2_old==2) = LL_old(Z2_old==2) - data.ACC_STEADY;
    LL_old(Z2_old==3) = LL_old(Z2_old==3) - data.ACC_EXPANSION;
    
    LL_old(Z1_old==1) = LL_old(Z1_old==1) - data.ACC_CONTRACTION;
    LL_old(Z1_old==2) = LL_old(Z1_old==2) - data.ACC_STEADY;
    LL_old(Z1_old==3) = LL_old(Z1_old==3) - data.ACC_EXPANSION;
    
    LL_new(Z2_new==1) = LL_new(Z2_new==1) - data.ACC_CONTRACTION;
    LL_new(Z2_new==2) = LL_new(Z2_new==2) - data.ACC_STEADY;
    LL_new(Z2_new==3) = LL_new(Z2_new==3) - data.ACC_EXPANSION;
    
    LL_new(Z1_new==1) = LL_new(Z1_new==1) - data.ACC_CONTRACTION;
    LL_new(Z1_new==2) = LL_new(Z1_new==2) - data.ACC_STEADY;
    LL_new(Z1_new==3) = LL_new(Z1_new==3) - data.ACC_EXPANSION;
    
    if ~strcmp(mode,'C14') && ~isnan(d18O(2,1))
        d18O_mu = interp1(target.stack(:,1),target.stack(:,2),A_old);
        d18O_stdv = interp1(target.stack(:,1),target.stack(:,3),A_old);
        NN = sum(~isnan(d18O(2,:)));
        for n = 1:NN
            % LL_d18O = (1-q)*exp(-(d18O(2,n)-d18O_mu).^2./(2*d18O_stdv.^2))./sqrt(2*pi*d18O_stdv.^2) + q/2*exp(-(d18O(2,n)-d18O_mu-d*d18O_stdv).^2./(2*d18O_stdv.^2))./sqrt(2*pi*d18O_stdv.^2) + q/2*exp(-(d18O(2,n)-d18O_mu+d*d18O_stdv).^2./(2*d18O_stdv.^2))./sqrt(2*pi*d18O_stdv.^2);
            % LL_old = LL_old + log(LL_d18O);
            LL_old = LL_old - (a+0.5)*log(1+(d18O(2,n)-d18O_mu).^2./(2*b*d18O_stdv.^2)) - log(d18O_stdv);
        end
        
        d18O_mu = interp1(target.stack(:,1),target.stack(:,2),A_new);
        d18O_stdv = interp1(target.stack(:,1),target.stack(:,3),A_new);
        for n = 1:NN
            % LL_d18O = (1-q)*exp(-(d18O(2,n)-d18O_mu).^2./(2*d18O_stdv.^2))./sqrt(2*pi*d18O_stdv.^2) + q/2*exp(-(d18O(2,n)-d18O_mu-d*d18O_stdv).^2./(2*d18O_stdv.^2))./sqrt(2*pi*d18O_stdv.^2) + q/2*exp(-(d18O(2,n)-d18O_mu+d*d18O_stdv).^2./(2*d18O_stdv.^2))./sqrt(2*pi*d18O_stdv.^2);
            % LL_new = LL_new + log(LL_d18O);
            LL_new = LL_new - (a+0.5)*log(1+(d18O(2,n)-d18O_mu).^2./(2*b*d18O_stdv.^2)) - log(d18O_stdv);
        end
    end
    
    if strcmp(mode,'d18O') == 0 && isempty(C14{2}) == 0
        for ll = 1:size(C14{2},1)
            c14_mu = interp1(target.cal_curve{C14{2}(ll,5)}(:,1),target.cal_curve{C14{2}(ll,5)}(:,2),A_old);
            c14_stdv = interp1(target.cal_curve{C14{2}(ll,5)}(:,1),target.cal_curve{C14{2}(ll,5)}(:,3),A_old);
            if C14{2}(ll,end) == 1
                LL_old = LL_old - 0.5*(c14_mu+C14{2}(ll,3)-C14{2}(ll,1)).^2./(c14_stdv.^2+C14{2}(ll,2)^2+C14{2}(ll,4)^2) - 0.5*log(c14_stdv.^2+C14{2}(ll,2)^2+C14{2}(ll,4)^2) - 0.5*log(2*pi);
            else
                LL_old = LL_old - (a_C14+0.5)*log(1+(c14_mu+C14{2}(ll,3)-C14{2}(ll,1)).^2./(2*b_C14*(c14_stdv.^2+C14{2}(ll,2)^2+C14{2}(ll,4)^2))) - 0.5*log(c14_stdv.^2+C14{2}(ll,2)^2+C14{2}(ll,4)^2);
            end
            
            c14_mu = interp1(target.cal_curve{C14{2}(ll,5)}(:,1),target.cal_curve{C14{2}(ll,5)}(:,2),A_new);
            c14_stdv = interp1(target.cal_curve{C14{2}(ll,5)}(:,1),target.cal_curve{C14{2}(ll,5)}(:,3),A_new);
            if C14{2}(ll,end) == 1
                LL_new = LL_new - 0.5*(c14_mu+C14{2}(ll,3)-C14{2}(ll,1)).^2./(c14_stdv.^2+C14{2}(ll,2)^2+C14{2}(ll,4)^2) - 0.5*log(c14_stdv.^2+C14{2}(ll,2)^2+C14{2}(ll,4)^2) - 0.5*log(2*pi);
            else
                LL_new = LL_new - (a_C14+0.5)*log(1+(c14_mu+C14{2}(ll,3)-C14{2}(ll,1)).^2./(2*b_C14*(c14_stdv.^2+C14{2}(ll,2)^2+C14{2}(ll,4)^2))) - 0.5*log(c14_stdv.^2+C14{2}(ll,2)^2+C14{2}(ll,4)^2);
            end
        end
    end
    
    if ~isnan(Age_Info(2,1))
        if Age_Info(2,3) == 0
            LL_old = LL_old - (A_old-Age_Info(2,1)).^2./(2*Age_Info(2,2).^2);
            LL_new = LL_new - (A_new-Age_Info(2,1)).^2./(2*Age_Info(2,2).^2);
        elseif Age_Info(2,3) == 1
            index_old = (A_old>Age_Info(2,1)-Age_Info(2,2))&(A_old<Age_Info(2,1)+Age_Info(2,2));
            index_new = (A_new>Age_Info(2,1)-Age_Info(2,2))&(A_new<Age_Info(2,1)+Age_Info(2,2));
            LL_old(~index_old) = -inf;
            LL_new(~index_new) = -inf;
        end
    end
    
    if ~isnan(data.max)
        LL_old(A_old>data.max) = -inf;
        LL_new(A_new>data.max) = -inf;
    end
    if ~isnan(data.min)
        LL_old(A_old<data.min) = -inf;
        LL_new(A_new<data.min) = -inf;
    end
    
    index = (LL_new-LL_old>Rand_Seed(2,:));
    A(2,index) = A_new(index);
    Z(1,index) = Z1_new(index);
    Z(2,index) = Z2_new(index);
    
    
    % rem(n,3) == 0:
    A_old = A(index_0,:);
    A_new = normrnd(A(index_0,:),(A(index_0_next,:)-A(index_0_prev,:))/8);
    
    Z_prev_old = Z(index_0_prev,:);
    Z_prev_new = 3*ones(N0,M);
    RR = interp1(R(:,1),R(:,2),A_new);
    Z_prev_new((A_new-A(index_0_prev,:))./(depth_diff(index_0_prev).*RR)<=1./0.9220) = 2;
    Z_prev_new((A_new-A(index_0_prev,:))./(depth_diff(index_0_prev).*RR)<=1./1.0850) = 1;
    
    Z_old = Z(index_0,:);
    Z_new = 3*ones(N0,M);
    RR = interp1(R(:,1),R(:,2),A(index_0_next,:));
    Z_new((A(index_0_next,:)-A_new)./(depth_diff(index_0).*RR)<=1./0.9220) = 2;
    Z_new((A(index_0_next,:)-A_new)./(depth_diff(index_0).*RR)<=1./1.0850) = 1;
    
    LL_old = zeros(N0,M);
    LL_new = zeros(N0,M);
    
    for k = 1:3
        for m = 1:3
            LL_old(Z(index_0_next,:)==k&Z_old==m) = PHI(k,m);
            LL_new(Z(index_0_next,:)==k&Z_new==m) = PHI(k,m);
        end
    end
    
    for k = 1:3
        for m = 1:3
            LL_old(Z_old==k&Z_prev_old==m) = LL_old(Z_old==k&Z_prev_old==m) + PHI(k,m);
            LL_new(Z_new==k&Z_prev_new==m) = LL_new(Z_new==k&Z_prev_new==m) + PHI(k,m);
        end
    end
    
    for k = 1:3
        for m = 1:3
            LL_old(Z_prev_old==k&Z(index_0_prev_prev,:)==m) = LL_old(Z_prev_old==k&Z(index_0_prev_prev,:)==m) + PHI(k,m);
            LL_new(Z_prev_new==k&Z(index_0_prev_prev,:)==m) = LL_new(Z_prev_new==k&Z(index_0_prev_prev,:)==m) + PHI(k,m);
        end
    end
    
    LL_new(A_new>=A(index_0_next,:)|A_new<=A(index_0_prev,:)) = -inf;
    
    TT_old = (A(index_0_next,:)-A_old)./(RR.*depth_diff(index_0));
    TT_new = (A(index_0_next,:)-A_new)./(RR.*depth_diff(index_0));
    
    LL_old(TT_old>0) = LL_old(TT_old>0) + interp1(data.ACC_MODEL(:,1),data.ACC_MODEL(:,2),TT_old(TT_old>0),'linear',-56);
    LL_new(TT_new>0) = LL_new(TT_new>0) + interp1(data.ACC_MODEL(:,1),data.ACC_MODEL(:,2),TT_new(TT_new>0),'linear',-56);
    
    LL_old(TT_old<=0) = -inf;
    LL_new(TT_new<=0) = -inf;
    
    LL_old(TT_old>1./data.lower_sedrate|TT_old<1./data.upper_sedrate) = -inf;
    LL_new(TT_new>1./data.lower_sedrate|TT_new<1./data.upper_sedrate) = -inf;
    
    RR_old = interp1(R(:,1),R(:,2),A_old);
    RR_new = interp1(R(:,1),R(:,2),A_new);
    
    TT_old = (A_old-A(index_0_prev,:))./(RR_old.*depth_diff(index_0_prev));
    TT_new = (A_new-A(index_0_prev,:))./(RR_new.*depth_diff(index_0_prev));
    
    LL_old(TT_old>0) = LL_old(TT_old>0) + interp1(data.ACC_MODEL(:,1),data.ACC_MODEL(:,2),TT_old(TT_old>0),'linear',-56);
    LL_new(TT_new>0) = LL_new(TT_new>0) + interp1(data.ACC_MODEL(:,1),data.ACC_MODEL(:,2),TT_new(TT_new>0),'linear',-56);
    
    LL_old(TT_old<=0) = -inf;
    LL_new(TT_new<=0) = -inf;
    
    LL_old(TT_old>1./data.lower_sedrate|TT_old<1./data.upper_sedrate) = -inf;
    LL_new(TT_new>1./data.lower_sedrate|TT_new<1./data.upper_sedrate) = -inf;
    
    LL_old(Z_old==1) = LL_old(Z_old==1) - data.ACC_CONTRACTION;
    LL_old(Z_old==2) = LL_old(Z_old==2) - data.ACC_STEADY;
    LL_old(Z_old==3) = LL_old(Z_old==3) - data.ACC_EXPANSION;
    
    LL_old(Z_prev_old==1) = LL_old(Z_prev_old==1) - data.ACC_CONTRACTION;
    LL_old(Z_prev_old==2) = LL_old(Z_prev_old==2) - data.ACC_STEADY;
    LL_old(Z_prev_old==3) = LL_old(Z_prev_old==3) - data.ACC_EXPANSION;
    
    LL_new(Z_new==1) = LL_new(Z_new==1) - data.ACC_CONTRACTION;
    LL_new(Z_new==2) = LL_new(Z_new==2) - data.ACC_STEADY;
    LL_new(Z_new==3) = LL_new(Z_new==3) - data.ACC_EXPANSION;
    
    LL_new(Z_prev_new==1) = LL_new(Z_prev_new==1) - data.ACC_CONTRACTION;
    LL_new(Z_prev_new==2) = LL_new(Z_prev_new==2) - data.ACC_STEADY;
    LL_new(Z_prev_new==3) = LL_new(Z_prev_new==3) - data.ACC_EXPANSION;
    
    if strcmp(mode,'C14') == 0
        d18O_mu = interp1(target.stack(:,1),target.stack(:,2),A_old);
        d18O_stdv = interp1(target.stack(:,1),target.stack(:,3),A_old);
        for n = 1:MAX
            % LL_d18O = (1-q)*exp(-(d18O(index_0,n)-d18O_mu).^2./(2*d18O_stdv.^2))./sqrt(2*pi*d18O_stdv.^2) + q/2*exp(-(d18O(index_0,n)-d18O_mu-d*d18O_stdv).^2./(2*d18O_stdv.^2))./sqrt(2*pi*d18O_stdv.^2) + q/2*exp(-(d18O(index_0,n)-d18O_mu+d*d18O_stdv).^2./(2*d18O_stdv.^2))./sqrt(2*pi*d18O_stdv.^2);
            % LL_d18O(isnan(LL_d18O)==1) = 1;
            % LL_old = LL_old + log(LL_d18O);
            LL_d18O =  - (a+0.5)*log(1+(d18O(index_0,n)-d18O_mu).^2./(2*b*d18O_stdv.^2)) - log(d18O_stdv);
            LL_d18O(isnan(LL_d18O)) = 0;
            LL_old = LL_old + LL_d18O;
        end
        
        d18O_mu = interp1(target.stack(:,1),target.stack(:,2),A_new);
        d18O_stdv = interp1(target.stack(:,1),target.stack(:,3),A_new);
        for n = 1:MAX
            % LL_d18O = (1-q)*exp(-(d18O(index_0,n)-d18O_mu).^2./(2*d18O_stdv.^2))./sqrt(2*pi*d18O_stdv.^2) + q/2*exp(-(d18O(index_0,n)-d18O_mu-d*d18O_stdv).^2./(2*d18O_stdv.^2))./sqrt(2*pi*d18O_stdv.^2) + q/2*exp(-(d18O(index_0,n)-d18O_mu+d*d18O_stdv).^2./(2*d18O_stdv.^2))./sqrt(2*pi*d18O_stdv.^2);
            % LL_d18O(isnan(LL_d18O)==1) = 1;
            % LL_new = LL_new + log(LL_d18O);
            LL_d18O =  - (a+0.5)*log(1+(d18O(index_0,n)-d18O_mu).^2./(2*b*d18O_stdv.^2)) - log(d18O_stdv);
            LL_d18O(isnan(LL_d18O)) = 0;
            LL_new = LL_new + LL_d18O;
        end
    end
    
    if strcmp(mode,'d18O') == 0
        LL_C14_old = zeros(N0,M);
        LL_C14_new = zeros(N0,M);
        m = 0;
        for n = 1:N
            if index_0(n) == 1
                m = m + 1;
                if isempty(C14{n}) == 0
                    for ll = 1:size(C14{n},1)
                        c14_mu = interp1(target.cal_curve{C14{n}(ll,5)}(:,1),target.cal_curve{C14{n}(ll,5)}(:,2),A_old(m,:));
                        c14_stdv = interp1(target.cal_curve{C14{n}(ll,5)}(:,1),target.cal_curve{C14{n}(ll,5)}(:,3),A_old(m,:));
                        if C14{n}(ll,end) == 1
                            LL_C14_old(m,:) = LL_C14_old(m,:) - 0.5*(c14_mu+C14{n}(ll,3)-C14{n}(ll,1)).^2./(c14_stdv.^2+C14{n}(ll,2)^2+C14{n}(ll,4)^2) - 0.5*log(c14_stdv.^2+C14{n}(ll,2)^2+C14{n}(ll,4)^2) - 0.5*log(2*pi);
                        else
                            LL_C14_old(m,:) = LL_C14_old(m,:) - (a_C14+0.5)*log(1+(c14_mu+C14{n}(ll,3)-C14{n}(ll,1)).^2./(2*b_C14*(c14_stdv.^2+C14{n}(ll,2)^2+C14{n}(ll,4)^2))) - 0.5*log(c14_stdv.^2+C14{n}(ll,2)^2+C14{n}(ll,4)^2);
                        end
                        
                        c14_mu = interp1(target.cal_curve{C14{n}(ll,5)}(:,1),target.cal_curve{C14{n}(ll,5)}(:,2),A_new(m,:));
                        c14_stdv = interp1(target.cal_curve{C14{n}(ll,5)}(:,1),target.cal_curve{C14{n}(ll,5)}(:,3),A_new(m,:));
                        if C14{n}(ll,end) == 1
                            LL_C14_new(m,:) = LL_C14_new(m,:) - 0.5*(c14_mu+C14{n}(ll,3)-C14{n}(ll,1)).^2./(c14_stdv.^2+C14{n}(ll,2)^2+C14{n}(ll,4)^2) - 0.5*log(c14_stdv.^2+C14{n}(ll,2)^2+C14{n}(ll,4)^2) - 0.5*log(2*pi);
                        else
                            LL_C14_new(m,:) = LL_C14_new(m,:) - (a_C14+0.5)*log(1+(c14_mu+C14{n}(ll,3)-C14{n}(ll,1)).^2./(2*b_C14*(c14_stdv.^2+C14{n}(ll,2)^2+C14{n}(ll,4)^2))) - 0.5*log(c14_stdv.^2+C14{n}(ll,2)^2+C14{n}(ll,4)^2);
                        end
                    end
                end
            end
        end
        LL_old = LL_old + LL_C14_old;
        LL_new = LL_new + LL_C14_new;
    end
    
    index_A = (Age_Info_0(:,3)==0);
    index_B = (Age_Info_0(:,3)==1);
    
    LL_old(index_A,:) = LL_old(index_A,:) - (A_old(index_A,:)-Age_Info_0(index_A,1)).^2./(2*Age_Info_0(index_A,2).^2);
    LL_new(index_A,:) = LL_new(index_A,:) - (A_new(index_A,:)-Age_Info_0(index_A,1)).^2./(2*Age_Info_0(index_A,2).^2);
    
    index_old = (A_old<=Age_Info_0(:,1)-Age_Info_0(:,2))|(A_old>=Age_Info_0(:,1)+Age_Info_0(:,2));
    index_new = (A_new<=Age_Info_0(:,1)-Age_Info_0(:,2))|(A_new>=Age_Info_0(:,1)+Age_Info_0(:,2));
    LL_old(index_old.*index_B==1) = -inf;
    LL_new(index_new.*index_B==1) = -inf;
    
    %{
    for n = 1:N0
        if isnan(Age_Info_0(n,1)) == 0
            if Age_Info_0(n,3) == 0
                LL_old(n,:) = LL_old(n,:) - (A_old(n,:)-Age_Info_0(n,1)).^2./(2*Age_Info_0(n,2).^2);
                LL_new(n,:) = LL_new(n,:) - (A_new(n,:)-Age_Info_0(n,1)).^2./(2*Age_Info_0(n,2).^2);
            elseif Age_Info_0(n,3) == 1
                index_old = (A_old(n,:)>Age_Info_0(n,1)-Age_Info_0(n,2))&(A_old(n,:)<Age_Info_0(n,1)+Age_Info_0(n,2));
                index_new = (A_new(n,:)>Age_Info_0(n,1)-Age_Info_0(n,2))&(A_new(n,:)<Age_Info_0(n,1)+Age_Info_0(n,2));
                LL_old(n,~index_old) = -inf;
                LL_new(n,~index_new) = -inf;
            end
        end
    end
    %}
    
    if ~isnan(data.max)
        LL_old(A_old>data.max) = -inf;
        LL_new(A_new>data.max) = -inf;
    end
    if ~isnan(data.min)
        LL_old(A_old<data.min) = -inf;
        LL_new(A_new<data.min) = -inf;
    end
    
    index = (LL_new-LL_old>Rand_Seed(index_0,:));
    AA = A(index_0,:);
    AA(index) = A_new(index);
    A(index_0,:) = AA;
    ZZ = Z(index_0,:);
    ZZ(index) = Z_new(index);
    Z(index_0,:) = ZZ;
    ZZ = Z(index_0_prev,:);
    ZZ(index) = Z_prev_new(index);
    Z(index_0_prev,:) = ZZ;
    
    
    % rem(n,3) == 1:
    A_old = A(index_1,:);
    A_new = normrnd(A(index_1,:),(A(index_1_next,:)-A(index_1_prev,:))/8);
    
    Z_prev_old = Z(index_1_prev,:);
    Z_prev_new = 3*ones(N1,M);
    RR = interp1(R(:,1),R(:,2),A_new);
    Z_prev_new((A_new-A(index_1_prev,:))./(depth_diff(index_1_prev).*RR)<=1./0.9220) = 2;
    Z_prev_new((A_new-A(index_1_prev,:))./(depth_diff(index_1_prev).*RR)<=1./1.0850) = 1;
    
    Z_old = Z(index_1,:);
    Z_new = 3*ones(N1,M);
    RR = interp1(R(:,1),R(:,2),A(index_1_next,:));
    Z_new((A(index_1_next,:)-A_new)./(depth_diff(index_1).*RR)<=1./0.9220) = 2;
    Z_new((A(index_1_next,:)-A_new)./(depth_diff(index_1).*RR)<=1./1.0850) = 1;
    
    LL_old = zeros(N1,M);
    LL_new = zeros(N1,M);
    
    for k = 1:3
        for m = 1:3
            LL_old(Z(index_1_next,:)==k&Z_old==m) = PHI(k,m);
            LL_new(Z(index_1_next,:)==k&Z_new==m) = PHI(k,m);
        end
    end
    
    for k = 1:3
        for m = 1:3
            LL_old(Z_old==k&Z_prev_old==m) = LL_old(Z_old==k&Z_prev_old==m) + PHI(k,m);
            LL_new(Z_new==k&Z_prev_new==m) = LL_new(Z_new==k&Z_prev_new==m) + PHI(k,m);
        end
    end
    
    for k = 1:3
        for m = 1:3
            LL_old(Z_prev_old==k&Z(index_1_prev_prev,:)==m) = LL_old(Z_prev_old==k&Z(index_1_prev_prev,:)==m) + PHI(k,m);
            LL_new(Z_prev_new==k&Z(index_1_prev_prev,:)==m) = LL_new(Z_prev_new==k&Z(index_1_prev_prev,:)==m) + PHI(k,m);
        end
    end
    
    LL_new(A_new>=A(index_1_next,:)|A_new<=A(index_1_prev,:)) = -inf;
    
    TT_old = (A(index_1_next,:)-A_old)./(RR.*depth_diff(index_1));
    TT_new = (A(index_1_next,:)-A_new)./(RR.*depth_diff(index_1));
    
    LL_old(TT_old>0) = LL_old(TT_old>0) + interp1(data.ACC_MODEL(:,1),data.ACC_MODEL(:,2),TT_old(TT_old>0),'linear',-56);
    LL_new(TT_new>0) = LL_new(TT_new>0) + interp1(data.ACC_MODEL(:,1),data.ACC_MODEL(:,2),TT_new(TT_new>0),'linear',-56);
    
    LL_old(TT_old<=0) = -inf;
    LL_new(TT_new<=0) = -inf;
    
    LL_old(TT_old>1./data.lower_sedrate|TT_old<1./data.upper_sedrate) = -inf;
    LL_new(TT_new>1./data.lower_sedrate|TT_new<1./data.upper_sedrate) = -inf;
    
    RR_old = interp1(R(:,1),R(:,2),A_old);
    RR_new = interp1(R(:,1),R(:,2),A_new);
    
    TT_old = (A_old-A(index_1_prev,:))./(RR_old.*depth_diff(index_1_prev));
    TT_new = (A_new-A(index_1_prev,:))./(RR_new.*depth_diff(index_1_prev));
    
    LL_old(TT_old>0) = LL_old(TT_old>0) + interp1(data.ACC_MODEL(:,1),data.ACC_MODEL(:,2),TT_old(TT_old>0),'linear',-56);
    LL_new(TT_new>0) = LL_new(TT_new>0) + interp1(data.ACC_MODEL(:,1),data.ACC_MODEL(:,2),TT_new(TT_new>0),'linear',-56);
    
    LL_old(TT_old<=0) = -inf;
    LL_new(TT_new<=0) = -inf;
    
    LL_old(TT_old>1./data.lower_sedrate|TT_old<1./data.upper_sedrate) = -inf;
    LL_new(TT_new>1./data.lower_sedrate|TT_new<1./data.upper_sedrate) = -inf;
    
    LL_old(Z_old==1) = LL_old(Z_old==1) - data.ACC_CONTRACTION;
    LL_old(Z_old==2) = LL_old(Z_old==2) - data.ACC_STEADY;
    LL_old(Z_old==3) = LL_old(Z_old==3) - data.ACC_EXPANSION;
    
    LL_old(Z_prev_old==1) = LL_old(Z_prev_old==1) - data.ACC_CONTRACTION;
    LL_old(Z_prev_old==2) = LL_old(Z_prev_old==2) - data.ACC_STEADY;
    LL_old(Z_prev_old==3) = LL_old(Z_prev_old==3) - data.ACC_EXPANSION;
    
    LL_new(Z_new==1) = LL_new(Z_new==1) - data.ACC_CONTRACTION;
    LL_new(Z_new==2) = LL_new(Z_new==2) - data.ACC_STEADY;
    LL_new(Z_new==3) = LL_new(Z_new==3) - data.ACC_EXPANSION;
    
    LL_new(Z_prev_new==1) = LL_new(Z_prev_new==1) - data.ACC_CONTRACTION;
    LL_new(Z_prev_new==2) = LL_new(Z_prev_new==2) - data.ACC_STEADY;
    LL_new(Z_prev_new==3) = LL_new(Z_prev_new==3) - data.ACC_EXPANSION;
    
    if strcmp(mode,'C14') == 0
        d18O_mu = interp1(target.stack(:,1),target.stack(:,2),A_old);
        d18O_stdv = interp1(target.stack(:,1),target.stack(:,3),A_old);
        for n = 1:MAX
            % LL_d18O = (1-q)*exp(-(d18O(index_1,n)-d18O_mu).^2./(2*d18O_stdv.^2))./sqrt(2*pi*d18O_stdv.^2) + q/2*exp(-(d18O(index_1,n)-d18O_mu-d*d18O_stdv).^2./(2*d18O_stdv.^2))./sqrt(2*pi*d18O_stdv.^2) + q/2*exp(-(d18O(index_1,n)-d18O_mu+d*d18O_stdv).^2./(2*d18O_stdv.^2))./sqrt(2*pi*d18O_stdv.^2);
            % LL_d18O(isnan(LL_d18O)==1) = 1;
            % LL_old = LL_old + log(LL_d18O);
            LL_d18O =  - (a+0.5)*log(1+(d18O(index_1,n)-d18O_mu).^2./(2*b*d18O_stdv.^2)) - log(d18O_stdv);
            LL_d18O(isnan(LL_d18O)) = 0;
            LL_old = LL_old + LL_d18O;
        end
        
        d18O_mu = interp1(target.stack(:,1),target.stack(:,2),A_new);
        d18O_stdv = interp1(target.stack(:,1),target.stack(:,3),A_new);
        for n = 1:MAX
            % LL_d18O = (1-q)*exp(-(d18O(index_1,n)-d18O_mu).^2./(2*d18O_stdv.^2))./sqrt(2*pi*d18O_stdv.^2) + q/2*exp(-(d18O(index_1,n)-d18O_mu-d*d18O_stdv).^2./(2*d18O_stdv.^2))./sqrt(2*pi*d18O_stdv.^2) + q/2*exp(-(d18O(index_1,n)-d18O_mu+d*d18O_stdv).^2./(2*d18O_stdv.^2))./sqrt(2*pi*d18O_stdv.^2);
            % LL_d18O(isnan(LL_d18O)==1) = 1;
            % LL_new = LL_new + log(LL_d18O);
            LL_d18O =  - (a+0.5)*log(1+(d18O(index_1,n)-d18O_mu).^2./(2*b*d18O_stdv.^2)) - log(d18O_stdv);
            LL_d18O(isnan(LL_d18O)) = 0;
            LL_new = LL_new + LL_d18O;
        end
    end
    
    if ~strcmp(mode,'d18O')
        LL_C14_old = zeros(N1,M);
        LL_C14_new = zeros(N1,M);
        m = 0;
        for n = 1:N
            if index_1(n) == 1
                m = m + 1;
                if isempty(C14{n}) == 0
                    for ll = 1:size(C14{n},1)
                        c14_mu = interp1(target.cal_curve{C14{n}(ll,5)}(:,1),target.cal_curve{C14{n}(ll,5)}(:,2),A_old(m,:));
                        c14_stdv = interp1(target.cal_curve{C14{n}(ll,5)}(:,1),target.cal_curve{C14{n}(ll,5)}(:,3),A_old(m,:));
                        if C14{n}(ll,end) == 1
                            LL_C14_old(m,:) = LL_C14_old(m,:) - 0.5*(c14_mu+C14{n}(ll,3)-C14{n}(ll,1)).^2./(c14_stdv.^2+C14{n}(ll,2)^2+C14{n}(ll,4)^2) - 0.5*log(c14_stdv.^2+C14{n}(ll,2)^2+C14{n}(ll,4)^2) - 0.5*log(2*pi);
                        else
                            LL_C14_old(m,:) = LL_C14_old(m,:) - (a_C14+0.5)*log(1+(c14_mu+C14{n}(ll,3)-C14{n}(ll,1)).^2./(2*b_C14*(c14_stdv.^2+C14{n}(ll,2)^2+C14{n}(ll,4)^2))) - 0.5*log(c14_stdv.^2+C14{n}(ll,2)^2+C14{n}(ll,4)^2);
                        end
                        
                        c14_mu = interp1(target.cal_curve{C14{n}(ll,5)}(:,1),target.cal_curve{C14{n}(ll,5)}(:,2),A_new(m,:));
                        c14_stdv = interp1(target.cal_curve{C14{n}(ll,5)}(:,1),target.cal_curve{C14{n}(ll,5)}(:,3),A_new(m,:));
                        if C14{n}(ll,end) == 1
                            LL_C14_new(m,:) = LL_C14_new(m,:) - 0.5*(c14_mu+C14{n}(ll,3)-C14{n}(ll,1)).^2./(c14_stdv.^2+C14{n}(ll,2)^2+C14{n}(ll,4)^2) - 0.5*log(c14_stdv.^2+C14{n}(ll,2)^2+C14{n}(ll,4)^2) - 0.5*log(2*pi);
                        else
                            LL_C14_new(m,:) = LL_C14_new(m,:) - (a_C14+0.5)*log(1+(c14_mu+C14{n}(ll,3)-C14{n}(ll,1)).^2./(2*b_C14*(c14_stdv.^2+C14{n}(ll,2)^2+C14{n}(ll,4)^2))) - 0.5*log(c14_stdv.^2+C14{n}(ll,2)^2+C14{n}(ll,4)^2);
                        end
                    end
                end
            end
        end
        LL_old = LL_old + LL_C14_old;
        LL_new = LL_new + LL_C14_new;
    end
    
    index_A = (Age_Info_1(:,3)==0);
    index_B = (Age_Info_1(:,3)==1);
    
    LL_old(index_A,:) = LL_old(index_A,:) - (A_old(index_A,:)-Age_Info_1(index_A,1)).^2./(2*Age_Info_1(index_A,2).^2);
    LL_new(index_A,:) = LL_new(index_A,:) - (A_new(index_A,:)-Age_Info_1(index_A,1)).^2./(2*Age_Info_1(index_A,2).^2);
    
    index_old = (A_old<=Age_Info_1(:,1)-Age_Info_1(:,2))|(A_old>=Age_Info_1(:,1)+Age_Info_1(:,2));
    index_new = (A_new<=Age_Info_1(:,1)-Age_Info_1(:,2))|(A_new>=Age_Info_1(:,1)+Age_Info_1(:,2));
    LL_old(index_old.*index_B==1) = -inf;
    LL_new(index_new.*index_B==1) = -inf;
    
    %{
    for n = 1:N1
        if isnan(Age_Info(index_1(n),1)) == 0
            if Age_Info(index_1(n),3) == 0
                LL_old(n,:) = LL_old(n,:) - (A_old(n,:)-Age_Info(index_1(n),1)).^2./(2*Age_Info(index_1(n),2).^2);
                LL_new(n,:) = LL_new(n,:) - (A_new(n,:)-Age_Info(index_1(n),1)).^2./(2*Age_Info(index_1(n),2).^2);
            elseif Age_Info(index_1(n),3) == 1
                index_old = (A_old(n,:)>Age_Info(index_1(n),1)-Age_Info(index_1(n),2))&(A_old(n,:)<Age_Info(index_1(n),1)+Age_Info(index_1(n),2));
                index_new = (A_new(n,:)>Age_Info(index_1(n),1)-Age_Info(index_1(n),2))&(A_new(n,:)<Age_Info(index_1(n),1)+Age_Info(index_1(n),2));
                LL_old(n,~index_old) = -inf;
                LL_new(n,~index_new) = -inf;
            end
        end
    end
    %}
    
    if ~isnan(data.max)
        LL_old(A_old>data.max) = -inf;
        LL_new(A_new>data.max) = -inf;
    end
    if ~isnan(data.min)
        LL_old(A_old<data.min) = -inf;
        LL_new(A_new<data.min) = -inf;
    end
    
    index = (LL_new-LL_old>Rand_Seed(index_1,:));
    AA = A(index_1,:);
    AA(index) = A_new(index);
    A(index_1,:) = AA;
    ZZ = Z(index_1,:);
    ZZ(index) = Z_new(index);
    Z(index_1,:) = ZZ;
    ZZ = Z(index_1_prev,:);
    ZZ(index) = Z_prev_new(index);
    Z(index_1_prev,:) = ZZ;
    
    
    % rem(n,3) == 2:
    A_old = A(index_2,:);
    A_new = normrnd(A(index_2,:),(A(index_2_next,:)-A(index_2_prev,:))/8);
    
    Z_prev_old = Z(index_2_prev,:);
    Z_prev_new = 3*ones(N2,M);
    RR = interp1(R(:,1),R(:,2),A_new);
    Z_prev_new((A_new-A(index_2_prev,:))./(depth_diff(index_2_prev).*RR)<=1./0.9220) = 2;
    Z_prev_new((A_new-A(index_2_prev,:))./(depth_diff(index_2_prev).*RR)<=1./1.0850) = 1;
    
    Z_old = Z(index_2,:);
    Z_new = 3*ones(N2,M);
    RR = interp1(R(:,1),R(:,2),A(index_2_next,:));
    Z_new((A(index_2_next,:)-A_new)./(depth_diff(index_2).*RR)<=1./0.9220) = 2;
    Z_new((A(index_2_next,:)-A_new)./(depth_diff(index_2).*RR)<=1./1.0850) = 1;
    
    LL_old = zeros(N2,M);
    LL_new = zeros(N2,M);
    
    for k = 1:3
        for m = 1:3
            LL_old(Z(index_2_next,:)==k&Z_old==m) = PHI(k,m);
            LL_new(Z(index_2_next,:)==k&Z_new==m) = PHI(k,m);
        end
    end
    
    for k = 1:3
        for m = 1:3
            LL_old(Z_old==k&Z_prev_old==m) = LL_old(Z_old==k&Z_prev_old==m) + PHI(k,m);
            LL_new(Z_new==k&Z_prev_new==m) = LL_new(Z_new==k&Z_prev_new==m) + PHI(k,m);
        end
    end
    
    for k = 1:3
        for m = 1:3
            LL_old(Z_prev_old==k&Z(index_2_prev_prev,:)==m) = LL_old(Z_prev_old==k&Z(index_2_prev_prev,:)==m) + PHI(k,m);
            LL_new(Z_prev_new==k&Z(index_2_prev_prev,:)==m) = LL_new(Z_prev_new==k&Z(index_2_prev_prev,:)==m) + PHI(k,m);
        end
    end
    
    LL_new(A_new>=A(index_2_next,:)|A_new<=A(index_2_prev,:)) = -inf;
    
    TT_old = (A(index_2_next,:)-A_old)./(RR.*depth_diff(index_2));
    TT_new = (A(index_2_next,:)-A_new)./(RR.*depth_diff(index_2));
    
    LL_old(TT_old>0) = LL_old(TT_old>0) + interp1(data.ACC_MODEL(:,1),data.ACC_MODEL(:,2),TT_old(TT_old>0),'linear',-56);
    LL_new(TT_new>0) = LL_new(TT_new>0) + interp1(data.ACC_MODEL(:,1),data.ACC_MODEL(:,2),TT_new(TT_new>0),'linear',-56);
    
    LL_old(TT_old<=0) = -inf;
    LL_new(TT_new<=0) = -inf;
    
    LL_old(TT_old>1./data.lower_sedrate|TT_old<1./data.upper_sedrate) = -inf;
    LL_new(TT_new>1./data.lower_sedrate|TT_new<1./data.upper_sedrate) = -inf;
    
    RR_old = interp1(R(:,1),R(:,2),A_old);
    RR_new = interp1(R(:,1),R(:,2),A_new);
    
    TT_old = (A_old-A(index_2_prev,:))./(RR_old.*depth_diff(index_2_prev));
    TT_new = (A_new-A(index_2_prev,:))./(RR_new.*depth_diff(index_2_prev));
    
    LL_old(TT_old>0) = LL_old(TT_old>0) + interp1(data.ACC_MODEL(:,1),data.ACC_MODEL(:,2),TT_old(TT_old>0),'linear',-56);
    LL_new(TT_new>0) = LL_new(TT_new>0) + interp1(data.ACC_MODEL(:,1),data.ACC_MODEL(:,2),TT_new(TT_new>0),'linear',-56);
    
    LL_old(TT_old<=0) = -inf;
    LL_new(TT_new<=0) = -inf;
    
    LL_old(TT_old>1./data.lower_sedrate|TT_old<1./data.upper_sedrate) = -inf;
    LL_new(TT_new>1./data.lower_sedrate|TT_new<1./data.upper_sedrate) = -inf;
    
    LL_old(Z_old==1) = LL_old(Z_old==1) - data.ACC_CONTRACTION;
    LL_old(Z_old==2) = LL_old(Z_old==2) - data.ACC_STEADY;
    LL_old(Z_old==3) = LL_old(Z_old==3) - data.ACC_EXPANSION;
    
    LL_old(Z_prev_old==1) = LL_old(Z_prev_old==1) - data.ACC_CONTRACTION;
    LL_old(Z_prev_old==2) = LL_old(Z_prev_old==2) - data.ACC_STEADY;
    LL_old(Z_prev_old==3) = LL_old(Z_prev_old==3) - data.ACC_EXPANSION;
    
    LL_new(Z_new==1) = LL_new(Z_new==1) - data.ACC_CONTRACTION;
    LL_new(Z_new==2) = LL_new(Z_new==2) - data.ACC_STEADY;
    LL_new(Z_new==3) = LL_new(Z_new==3) - data.ACC_EXPANSION;
    
    LL_new(Z_prev_new==1) = LL_new(Z_prev_new==1) - data.ACC_CONTRACTION;
    LL_new(Z_prev_new==2) = LL_new(Z_prev_new==2) - data.ACC_STEADY;
    LL_new(Z_prev_new==3) = LL_new(Z_prev_new==3) - data.ACC_EXPANSION;
    
    if strcmp(mode,'C14') == 0
        d18O_mu = interp1(target.stack(:,1),target.stack(:,2),A_old);
        d18O_stdv = interp1(target.stack(:,1),target.stack(:,3),A_old);
        for n = 1:MAX
            % LL_d18O = (1-q)*exp(-(d18O(index_2,n)-d18O_mu).^2./(2*d18O_stdv.^2))./sqrt(2*pi*d18O_stdv.^2) + q/2*exp(-(d18O(index_2,n)-d18O_mu-d*d18O_stdv).^2./(2*d18O_stdv.^2))./sqrt(2*pi*d18O_stdv.^2) + q/2*exp(-(d18O(index_2,n)-d18O_mu+d*d18O_stdv).^2./(2*d18O_stdv.^2))./sqrt(2*pi*d18O_stdv.^2);
            % LL_d18O(isnan(LL_d18O)==1) = 1;
            % LL_old = LL_old + log(LL_d18O);
            LL_d18O =  - (a+0.5)*log(1+(d18O(index_2,n)-d18O_mu).^2./(2*b*d18O_stdv.^2)) - log(d18O_stdv);
            LL_d18O(isnan(LL_d18O)) = 0;
            LL_old = LL_old + LL_d18O;
        end
        
        d18O_mu = interp1(target.stack(:,1),target.stack(:,2),A_new);
        d18O_stdv = interp1(target.stack(:,1),target.stack(:,3),A_new);
        for n = 1:MAX
            % LL_d18O = (1-q)*exp(-(d18O(index_2,n)-d18O_mu).^2./(2*d18O_stdv.^2))./sqrt(2*pi*d18O_stdv.^2) + q/2*exp(-(d18O(index_2,n)-d18O_mu-d*d18O_stdv).^2./(2*d18O_stdv.^2))./sqrt(2*pi*d18O_stdv.^2) + q/2*exp(-(d18O(index_2,n)-d18O_mu+d*d18O_stdv).^2./(2*d18O_stdv.^2))./sqrt(2*pi*d18O_stdv.^2);
            % LL_d18O(isnan(LL_d18O)==1) = 1;
            % LL_new = LL_new + log(LL_d18O);
            LL_d18O =  - (a+0.5)*log(1+(d18O(index_2,n)-d18O_mu).^2./(2*b*d18O_stdv.^2)) - log(d18O_stdv);
            LL_d18O(isnan(LL_d18O)) = 0;
            LL_new = LL_new + LL_d18O;
        end
    end
    
    if strcmp(mode,'d18O') == 0
        LL_C14_old = zeros(N2,M);
        LL_C14_new = zeros(N2,M);
        m = 0;
        for n = 1:N
            if index_2(n) == 1
                m = m + 1;
                if isempty(C14{n}) == 0
                    for ll = 1:size(C14{n},1)
                        c14_mu = interp1(target.cal_curve{C14{n}(ll,5)}(:,1),target.cal_curve{C14{n}(ll,5)}(:,2),A_old(m,:));
                        c14_stdv = interp1(target.cal_curve{C14{n}(ll,5)}(:,1),target.cal_curve{C14{n}(ll,5)}(:,3),A_old(m,:));
                        if C14{n}(ll,end) == 1
                            LL_C14_old(m,:) = LL_C14_old(m,:) - 0.5*(c14_mu+C14{n}(ll,3)-C14{n}(ll,1)).^2./(c14_stdv.^2+C14{n}(ll,2)^2+C14{n}(ll,4)^2) - 0.5*log(c14_stdv.^2+C14{n}(ll,2)^2+C14{n}(ll,4)^2) - 0.5*log(2*pi);
                        else
                            LL_C14_old(m,:) = LL_C14_old(m,:) - (a_C14+0.5)*log(1+(c14_mu+C14{n}(ll,3)-C14{n}(ll,1)).^2./(2*b_C14*(c14_stdv.^2+C14{n}(ll,2)^2+C14{n}(ll,4)^2))) - 0.5*log(c14_stdv.^2+C14{n}(ll,2)^2+C14{n}(ll,4)^2);
                        end
                        
                        c14_mu = interp1(target.cal_curve{C14{n}(ll,5)}(:,1),target.cal_curve{C14{n}(ll,5)}(:,2),A_new(m,:));
                        c14_stdv = interp1(target.cal_curve{C14{n}(ll,5)}(:,1),target.cal_curve{C14{n}(ll,5)}(:,3),A_new(m,:));
                        if C14{n}(ll,end) == 1
                            LL_C14_new(m,:) = LL_C14_new(m,:) - 0.5*(c14_mu+C14{n}(ll,3)-C14{n}(ll,1)).^2./(c14_stdv.^2+C14{n}(ll,2)^2+C14{n}(ll,4)^2) - 0.5*log(c14_stdv.^2+C14{n}(ll,2)^2+C14{n}(ll,4)^2) - 0.5*log(2*pi);
                        else
                            LL_C14_new(m,:) = LL_C14_new(m,:) - (a_C14+0.5)*log(1+(c14_mu+C14{n}(ll,3)-C14{n}(ll,1)).^2./(2*b_C14*(c14_stdv.^2+C14{n}(ll,2)^2+C14{n}(ll,4)^2))) - 0.5*log(c14_stdv.^2+C14{n}(ll,2)^2+C14{n}(ll,4)^2);
                        end
                    end
                end
            end
        end
        LL_old = LL_old + LL_C14_old;
        LL_new = LL_new + LL_C14_new;
    end
    
    index_A = (Age_Info_2(:,3)==0);
    index_B = (Age_Info_2(:,3)==1);
    
    LL_old(index_A,:) = LL_old(index_A,:) - (A_old(index_A,:)-Age_Info_2(index_A,1)).^2./(2*Age_Info_2(index_A,2).^2);
    LL_new(index_A,:) = LL_new(index_A,:) - (A_new(index_A,:)-Age_Info_2(index_A,1)).^2./(2*Age_Info_2(index_A,2).^2);
    
    index_old = (A_old<=Age_Info_2(:,1)-Age_Info_2(:,2))|(A_old>=Age_Info_2(:,1)+Age_Info_2(:,2));
    index_new = (A_new<=Age_Info_2(:,1)-Age_Info_2(:,2))|(A_new>=Age_Info_2(:,1)+Age_Info_2(:,2));
    LL_old(index_old.*index_B==1) = -inf;
    LL_new(index_new.*index_B==1) = -inf;
    
    %{
    for n = 1:N2
        if isnan(Age_Info(index_2(n),1)) == 0
            if Age_Info(index_2(n),3) == 0
                LL_old(n,:) = LL_old(n,:) - (A_old(n,:)-Age_Info(index_2(n),1)).^2./(2*Age_Info(index_2(n),2).^2);
                LL_new(n,:) = LL_new(n,:) - (A_new(n,:)-Age_Info(index_2(n),1)).^2./(2*Age_Info(index_2(n),2).^2);
            elseif Age_Info(index_2(n),3) == 1
                index_old = (A_old(n,:)>Age_Info(index_2(n),1)-Age_Info(index_2(n),2))&(A_old(n,:)<Age_Info(index_2(n),1)+Age_Info(index_2(n),2));
                index_new = (A_new(n,:)>Age_Info(index_2(n),1)-Age_Info(index_2(n),2))&(A_new(n,:)<Age_Info(index_2(n),1)+Age_Info(index_2(n),2));
                LL_old(n,~index_old) = -inf;
                LL_new(n,~index_new) = -inf;
            end
        end
    end
    %}
    
    if ~isnan(data.max)
        LL_old(A_old>data.max) = -inf;
        LL_new(A_new>data.max) = -inf;
    end
    if ~isnan(data.min)
        LL_old(A_old<data.min) = -inf;
        LL_new(A_new<data.min) = -inf;
    end
    
    index = (LL_new-LL_old>Rand_Seed(index_2,:));
    AA = A(index_2,:);
    AA(index) = A_new(index);
    A(index_2,:) = AA;
    ZZ = Z(index_2,:);
    ZZ(index) = Z_new(index);
    Z(index_2,:) = ZZ;
    ZZ = Z(index_2_prev,:);
    ZZ(index) = Z_prev_new(index);
    Z(index_2_prev,:) = ZZ;
    
    
    % n == N:
    A_old = A(N,:);
    A_new = normrnd(A(N,:),(A(N,:)-A(N-1,:)+1)/8);
    
    ZN1_old = Z(N-1,:);
    ZN1_new = 3*ones(1,M);
    RR = interp1(R(:,1),R(:,2),A_new);
    ZN1_new((A_new-A(N-1,:))./(depth_diff(N-1).*RR)<=1./0.9220) = 2;
    ZN1_new((A_new-A(N-1,:))./(depth_diff(N-1).*RR)<=1./1.0850) = 1;
    
    LL_old = zeros(1,M);
    LL_new = zeros(1,M);
    
    for m = 1:3
        LL_old(ZN1_old==m) = PHI(4,m);
        LL_new(ZN1_new==m) = PHI(4,m);
    end
    
    for k = 1:3
        for m = 1:3
            LL_old(ZN1_old==k&Z(N-2,:)==m) = PHI(k,m);
            LL_new(ZN1_new==k&Z(N-2,:)==m) = PHI(k,m);
        end
    end
    
    LL_new(A_new<=A(N-1,:)|A_new>target.stack(end,1)) = -inf;
    
    RR_old = interp1(R(:,1),R(:,2),A_old);
    RR_new = interp1(R(:,1),R(:,2),A_new);
    TT_old = (A_old-A(N-1,:))./(RR_old.*depth_diff(N-1));
    TT_new = (A_new-A(N-1,:))./(RR_new.*depth_diff(N-1));
    
    LL_old(TT_old>0) = LL_old(TT_old>0) + interp1(data.ACC_MODEL(:,1),data.ACC_MODEL(:,2),TT_old(TT_old>0),'linear',-56);
    LL_new(TT_new>0) = LL_new(TT_new>0) + interp1(data.ACC_MODEL(:,1),data.ACC_MODEL(:,2),TT_new(TT_new>0),'linear',-56);
    
    LL_old(TT_old<=0) = -inf;
    LL_new(TT_new<=0) = -inf;
    
    LL_old(TT_old>1./data.lower_sedrate|TT_old<1./data.upper_sedrate) = -inf;
    LL_new(TT_new>1./data.lower_sedrate|TT_new<1./data.upper_sedrate) = -inf;
    
    LL_new = LL_new - 0.5*(A_new-A_old).^2./((A_new-A(N-1,:)+1).^2/64) - log((A_new-A(N-1,:)+1)/8);
    LL_old = LL_old - 0.5*(A_new-A_old).^2./((A_old-A(N-1,:)+1).^2/64) - log((A_old-A(N-1,:)+1)/8);
    
    LL_old(ZN1_old==1) = LL_old(ZN1_old==1) - data.ACC_CONTRACTION;
    LL_old(ZN1_old==2) = LL_old(ZN1_old==2) - data.ACC_STEADY;
    LL_old(ZN1_old==3) = LL_old(ZN1_old==3) - data.ACC_EXPANSION;
    
    LL_new(ZN1_new==1) = LL_new(ZN1_new==1) - data.ACC_CONTRACTION;
    LL_new(ZN1_new==2) = LL_new(ZN1_new==2) - data.ACC_STEADY;
    LL_new(ZN1_new==3) = LL_new(ZN1_new==3) - data.ACC_EXPANSION;
    
    if ~strcmp(mode,'C14') && ~isnan(d18O(N,1))
        d18O_mu = interp1(target.stack(:,1),target.stack(:,2),A_old);
        d18O_stdv = interp1(target.stack(:,1),target.stack(:,3),A_old);
        NN = sum(~isnan(d18O(N,:)));
        for n = 1:NN
            % LL_d18O = (1-q)*exp(-(d18O(N,n)-d18O_mu).^2./(2*d18O_stdv.^2))./sqrt(2*pi*d18O_stdv.^2) + q/2*exp(-(d18O(N,n)-d18O_mu-d*d18O_stdv).^2./(2*d18O_stdv.^2))./sqrt(2*pi*d18O_stdv.^2) + q/2*exp(-(d18O(N,n)-d18O_mu+d*d18O_stdv).^2./(2*d18O_stdv.^2))./sqrt(2*pi*d18O_stdv.^2);
            % LL_old = LL_old + log(LL_d18O);
            LL_old = LL_old - (a+0.5)*log(1+(d18O(N,n)-d18O_mu).^2./(2*b*d18O_stdv.^2)) - log(d18O_stdv);
        end
        
        d18O_mu = interp1(target.stack(:,1),target.stack(:,2),A_new);
        d18O_stdv = interp1(target.stack(:,1),target.stack(:,3),A_new);
        for n = 1:NN
            % LL_d18O = (1-q)*exp(-(d18O(N,n)-d18O_mu).^2./(2*d18O_stdv.^2))./sqrt(2*pi*d18O_stdv.^2) + q/2*exp(-(d18O(N,n)-d18O_mu-d*d18O_stdv).^2./(2*d18O_stdv.^2))./sqrt(2*pi*d18O_stdv.^2) + q/2*exp(-(d18O(N,n)-d18O_mu+d*d18O_stdv).^2./(2*d18O_stdv.^2))./sqrt(2*pi*d18O_stdv.^2);
            % LL_new = LL_new + log(LL_d18O);
            LL_new = LL_new - (a+0.5)*log(1+(d18O(N,n)-d18O_mu).^2./(2*b*d18O_stdv.^2)) - log(d18O_stdv);
        end
    end
    
    if ~strcmp(mode,'d18O') && ~isempty(C14{N})
        for ll = 1:size(C14{N},1)
            c14_mu = interp1(target.cal_curve{C14{N}(ll,5)}(:,1),target.cal_curve{C14{N}(ll,5)}(:,2),A_old);
            c14_stdv = interp1(target.cal_curve{C14{N}(ll,5)}(:,1),target.cal_curve{C14{N}(ll,5)}(:,3),A_old);
            if C14{N}(ll,end) == 1
                LL_old = LL_old - 0.5*(c14_mu+C14{N}(ll,3)-C14{N}(ll,1)).^2./(c14_stdv.^2+C14{N}(ll,2)^2+C14{N}(ll,4)^2) - 0.5*log(c14_stdv.^2+C14{N}(ll,2)^2+C14{N}(ll,4)^2) - 0.5*log(2*pi);
            else
                LL_old = LL_old - (a_C14+0.5)*log(1+(c14_mu+C14{N}(ll,3)-C14{N}(ll,1)).^2./(2*b_C14*(c14_stdv.^2+C14{N}(ll,2)^2+C14{N}(ll,4)^2))) - 0.5*log(c14_stdv.^2+C14{N}(ll,2)^2+C14{N}(ll,4)^2);
            end
            
            c14_mu = interp1(target.cal_curve{C14{N}(ll,5)}(:,1),target.cal_curve{C14{N}(ll,5)}(:,2),A_new);
            c14_stdv = interp1(target.cal_curve{C14{N}(ll,5)}(:,1),target.cal_curve{C14{N}(ll,5)}(:,3),A_new);
            if C14{N}(ll,end) == 1
                LL_new = LL_new - 0.5*(c14_mu+C14{N}(ll,3)-C14{N}(ll,1)).^2./(c14_stdv.^2+C14{N}(ll,2)^2+C14{N}(ll,4)^2) - 0.5*log(c14_stdv.^2+C14{N}(ll,2)^2+C14{N}(ll,4)^2) - 0.5*log(2*pi);
            else
                LL_new = LL_new - (a_C14+0.5)*log(1+(c14_mu+C14{N}(ll,3)-C14{N}(ll,1)).^2./(2*b_C14*(c14_stdv.^2+C14{N}(ll,2)^2+C14{N}(ll,4)^2))) - 0.5*log(c14_stdv.^2+C14{N}(ll,2)^2+C14{N}(ll,4)^2);
            end
        end
    end
    
    if ~isnan(Age_Info(N,1))
        if Age_Info(N,3) == 0
            LL_old = LL_old - (A_old-Age_Info(N,1)).^2./(2*Age_Info(N,2).^2);
            LL_new = LL_new - (A_new-Age_Info(N,1)).^2./(2*Age_Info(N,2).^2);
        elseif Age_Info(N,3) == 1
            index_old = (A_old>Age_Info(N,1)-Age_Info(N,2))&(A_old<Age_Info(N,1)+Age_Info(N,2));
            index_new = (A_new>Age_Info(N,1)-Age_Info(N,2))&(A_new<Age_Info(N,1)+Age_Info(N,2));
            LL_old(~index_old) = -inf;
            LL_new(~index_new) = -inf;
        end
    end
    
    if ~isnan(data.max)
        LL_old(A_old>data.max) = -inf;
        LL_new(A_new>data.max) = -inf;
    end
    if ~isnan(data.min)
        LL_old(A_old<data.min) = -inf;
        LL_new(A_new<data.min) = -inf;
    end
    
    index = (LL_new-LL_old>Rand_Seed(N,:));
    A(N,index) = A_new(index);
    Z(N-1,index) = ZN1_new(index);
end


HH = cell(M,1);
if strcmp(mode,'C14') == 0
    for m = 1:M
        rand_seed = rand(N,MAX);
        
        mu = interp1(target.stack(:,1),target.stack(:,2),A(:,m));
        sig = interp1(target.stack(:,1),target.stack(:,3),A(:,m));
        
        DET0 = log(1-q) - (d18O-mu).^2./(2*sig.^2) - log(sig);
        DET1 = log(q/2) + log(exp(-(d18O-mu-d*sig).^2./(2*sig.^2))+exp(-(d18O-mu+d*sig).^2./(2*sig.^2))) - log(sig);
        
        AMAX = max(DET0,DET1);
        DET0 = exp(DET0-AMAX);
        DET1 = exp(DET1-AMAX);
        DET = DET1./(DET0+DET1);
        HH{m} = (rand_seed<DET).*1;
        HH{m}(isnan(d18O)) = 1;
    end
end


Samples.ages = A;
Samples.isoutlier = HH;


QQ = zeros(N,S);
index = ceil(M*rand(N,S));
for n = 1:N
    QQ(n,:) = A(n,index(n,:));
end


end