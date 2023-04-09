function [param,target] = updateStack(data,Samples,param,target,setting)

L = length(Samples);
T = size(target.stack,1);

X_min = target.stack(1,1);
X_max = target.stack(end,1);

AGE = target.stack(:,1);

M = size(Samples(1).ages,2);

X0 = cell(M,1);

for m = 1:M
    X0{m} = (floor(X_min):setting.interval_induced:ceil(X_max))';
    rand_seed = (2*rand(1,M)-1)*setting.interval_induced;
    X0{m} = X0{m} + rand_seed;
end



X = cell(M,1);
Y = cell(M,1);
H = cell(M,1);
R = cell(M,1);

for m = 1:M
    XX = cell(L,1);
    YY = cell(L,1);
    HH = cell(L,1);
    
    for ll = 1:L
        K = size(data(ll).d18O,2);
        AA = repmat(Samples(ll).ages(:,m),[1,K]);
        BB = (data(ll).d18O-data(ll).shift-data(ll).scale*param.shift)/(data(ll).scale*param.scale);
        CC = Samples(ll).isoutlier{m};
        
        index = ~isnan(BB);
        AA = AA(index);
        BB = BB(index);
        CC = CC(index);
        
        XX{ll} = AA;
        YY{ll} = BB;
        HH{ll} = CC;
    end
    
    X{m} = vertcat(XX{:});
    Y{m} = vertcat(YY{:});
    H{m} = vertcat(HH{:});
end

BIAS = (quantile(Y{m},0.9)+quantile(Y{m},0.1))/2;

for m = 1:M
    IND = zeros(size(X0{m},1),1);
    for k = 1:length(IND)
        if sum(abs(X0{m}(k)-X{m})<setting.interval_induced/2) > 0
            IND(k) = 1;
        end
    end
    X0{m} = X0{m}(IND==1);
end

for m = 1:M
    N = size(X{m},1);
    R{m} = param.lambda*ones(N,1);
    Y{m} = Y{m} - BIAS;
end


% Learning parameters:
for r = 1:10
    disp(['   Iteration ',num2str(r),'/10...']);
    % Constructing regression models:
    [MU,NU] = learnParam_Reg(X0,X,Y,H,R,setting.kernel_function,param);
    % Learning variances:
    [param,R] = learnParam_Var(X,Y,H,MU,NU,param,setting.variance,r);
    % Sampling outliers:
    H = learnParam_Outlier(X,Y,MU,NU,R,param);
    % Learning parameters:
    param = learnParam(X0,X,Y,H,R,setting.kernel_function,param);
end
[MU,NU] = learnParam_Reg(X0,X,Y,H,R,setting.kernel_function,param);
for m = 1:M
    param.mu{m} = MU{m};
    param.nu{m} = NU{m};
    Y{m} = Y{m} + BIAS;
end



MU = zeros(T,M);
NU = zeros(T,M);

STACK_SAMPLE = zeros(T,M,2);


TT = AGE;

K_TT = diag(getCov(TT,TT,param,setting.kernel_function)) + param.lambda0^2;

for m = 1:M
    N0 = size(X0{m},1);
    
    XX0 = X0{m};
    XX = X{m}(H{m}==0);
    YY = Y{m}(H{m}==0) - BIAS;
    
    N = size(XX,1);
    if strcmp(setting.variance,'heteroscedastic') == 1
        QQ = abs(XX-XX');
        PP = sort(QQ,2,'ascend');
        BW = PP(:,ceil(N*param.K(m)/100));
        
        QQ = abs(AGE-XX');
        
        WW = - 0.5*QQ.^2./(BW'.^2) - log(BW');
        WW = WW - max(WW,[],2);
        WW = exp(WW);
        LL = WW./sum(WW,2);
        
        RRT = sum(((YY'-param.mu{m}(H{m}==0)').^2+param.nu{m}(H{m}==0)').*LL,2);
    elseif strcmp(setting.variance,'homoscedastic') == 1
        RRT = param.lambda^2*ones(T,1);
    end
    
    RR = interp1(AGE,RRT,XX);
    
    K_00 = getCov(XX0,XX0,param,setting.kernel_function) + param.lambda0^2*eye(N0);
    K_X0 = getCov(XX,XX0,param,setting.kernel_function);
    K_T0 = getCov(TT,XX0,param,setting.kernel_function);
    
    WW_inv = (K_00+K_X0'*(K_X0./RR))\eye(N0);
    QQ_inv = WW_inv - K_00\eye(N0);
    
    MU(:,m) = K_T0*(WW_inv*(K_X0'*(YY./RR))) + BIAS;
    NU(:,m) = K_TT + RRT + sum((K_T0*QQ_inv).*K_T0,2);
    
    % MU(:,m) = MU(:,m)*param.scale + param.shift;
    % NU(:,m) = NU(:,m)*param.scale^2;
    
    SS = getCov(TT,TT,param,setting.kernel_function) + param.lambda0^2*eye(T);
    SS = SS + K_T0*QQ_inv*K_T0';
    SS = (SS+SS')/2;
    % STACK_SAMPLE(:,m,1) = MU(:,m) + param.scale*mvnrnd(zeros(1,T),SS)';
    STACK_SAMPLE(:,m,1) = MU(:,m) + mvnrnd(zeros(1,T),SS)';
    
    SS = SS + diag(RRT);
    SS = (SS+SS')/2;
    % STACK_SAMPLE(:,m,2) = MU(:,m) + param.scale*mvnrnd(zeros(1,T),SS)';
    STACK_SAMPLE(:,m,2) = MU(:,m) + mvnrnd(zeros(1,T),SS)';
end


target.stack(:,2) = mean(MU,2);
target.stack(:,3) = sqrt(mean(NU+(MU-target.stack(:,2)).^2,2));

target.stack_sample = STACK_SAMPLE;


end