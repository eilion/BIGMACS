function [param] = learnParam(X0,X,Y,H,R,kernel,param)

beta1 = 0.9;
beta2 = 0.999;
epsilon = 1e-8;
gamma = 1e-3;
eps = 1e-4;
max_iters = 3000;

M = length(X);

SCH = ceil(M*rand(1,max_iters));
SCH2 = zeros(M,max_iters);

for m = 1:M
    NN0 = size(X0{m},1);
    if NN0 > 100
        QQ = 100;
        WW = median(X0{m},2);
        VV = sort(X{m},'ascend');
        for n = 101:NN0
            if sum(VV>WW(n-100)) > 100
                QQ = n;
            end
        end
        QQ = min(QQ);
        SCH2(m,:) = ceil((QQ-99)*rand(1,max_iters));
    end
end

for m = 1:M
    X{m} = X{m}(H{m}==0,:);
    Y{m} = Y{m}(H{m}==0,:);
    R{m} = R{m}(H{m}==0,:);
end


Mw = zeros(3,1);
Vw = zeros(3,1);

for r = 1:max_iters
    
    PDEV = zeros(3,1);
    
    m = SCH(r);
    NN0 = size(X0{m},1);
    
    if NN0 > 100
        XX0 = X0{m}(SCH2(m,r):SCH2(m,r)+99);
        N0 = 100;
    else
        XX0 = X0{m};
        N0 = NN0;
    end
    XX = X{m};
    YY = Y{m};
    RR = R{m};
    
    index = ((XX>XX0(1))&(XX<XX0(end)));
    XX = XX(index,:);
    YY = YY(index,:);
    RR = RR(index,:);
    
    N = size(XX,1);
    
    % eta:
    RR0 = param.lambda^2.*RR;
    
    K_00 = getCov(XX0,XX0,param,kernel) + param.lambda0^2*eye(N0);
    K_00_inv = K_00\eye(N0);
    K_X0 = getCov(XX,XX0,param,kernel);
    
    WW_inv = (K_00+K_X0'*(K_X0./RR0))\eye(N0);
    GG = YY./RR0;
    
    BB = - 0.5*YY'*GG + 0.5*(K_X0'*GG)'*WW_inv*(K_X0'*GG);
    AA = param.eta^2 + param.lambda0^2 - sum((K_X0*K_00_inv).*K_X0,2);
    
    LL0 = BB - 0.5*log(det(eye(N0)+K_X0'*(K_X0./RR0)*K_00_inv)) - 0.5*sum(AA./RR0) - N*log(param.lambda);
    
    param_new = param;
    param_new.eta = param.eta + eps;
    
    K_00 = getCov(XX0,XX0,param_new,kernel) + param_new.lambda0^2*eye(N0);
    K_00_inv = K_00\eye(N0);
    K_X0 = getCov(XX,XX0,param_new,kernel);
    
    WW_inv = (K_00+K_X0'*(K_X0./RR0))\eye(N0);
    GG = YY./RR0;
    
    BB = - 0.5*YY'*GG + 0.5*(K_X0'*GG)'*WW_inv*(K_X0'*GG);
    AA = param_new.eta^2 + param_new.lambda0^2 - sum((K_X0*K_00_inv).*K_X0,2);
    
    LL1 = BB - 0.5*log(det(eye(N0)+K_X0'*(K_X0./RR0)*K_00_inv)) - 0.5*sum(AA./RR0) - N*log(param_new.lambda);
    
    PDEV(1) = (LL1-LL0)/eps;
    
    % xi:
    param_new = param;
    param_new.xi = param.xi + eps;
    
    K_00 = getCov(XX0,XX0,param_new,kernel) + param_new.lambda0^2*eye(N0);
    K_00_inv = K_00\eye(N0);
    K_X0 = getCov(XX,XX0,param_new,kernel);
    
    WW_inv = (K_00+K_X0'*(K_X0./RR0))\eye(N0);
    GG = YY./RR0;
    
    BB = - 0.5*YY'*GG + 0.5*(K_X0'*GG)'*WW_inv*(K_X0'*GG);
    AA = param_new.eta^2 + param_new.lambda0^2 - sum((K_X0*K_00_inv).*K_X0,2);
    
    LL1 = BB - 0.5*log(det(eye(N0)+K_X0'*(K_X0./RR0)*K_00_inv)) - 0.5*sum(AA./RR0) - N*log(param_new.lambda);
    
    PDEV(2) = (LL1-LL0)/eps;
    
    % if strcmp(VAR,'homoscedastic') == 1
    % lambda:
    param_new = param;
    param_new.lambda = param.lambda + eps;
    
    RR1 = param_new.lambda^2.*RR;
    
    K_00 = getCov(XX0,XX0,param_new,kernel) + param_new.lambda0^2*eye(N0);
    K_00_inv = K_00\eye(N0);
    K_X0 = getCov(XX,XX0,param_new,kernel);
    
    WW_inv = (K_00+K_X0'*(K_X0./RR1))\eye(N0);
    GG = YY./RR1;
    
    BB = - 0.5*YY'*GG + 0.5*(K_X0'*GG)'*WW_inv*(K_X0'*GG);
    AA = param_new.eta^2 + param_new.lambda0^2 - sum((K_X0*K_00_inv).*K_X0,2);
    
    LL1 = BB - 0.5*log(det(eye(N0)+K_X0'*(K_X0./RR1)*K_00_inv)) - 0.5*sum(AA./RR1) - N*log(param_new.lambda);
    
    PDEV(3) = (LL1-LL0)/eps;
    % end
    
    
    
    Mw = beta1*Mw - (1-beta1)*PDEV;
    Vw = beta2*Vw + (1-beta2)*PDEV.*PDEV;
    
    param.eta = param.eta - (gamma*sqrt(1-beta2^r)/(1-beta1^r)).*Mw(1)./(sqrt(Vw(1))+epsilon);
    param.xi = param.xi - (gamma*sqrt(1-beta2^r)/(1-beta1^r)).*Mw(2)./(sqrt(Vw(2))+epsilon);
    param.lambda = param.lambda - (gamma*sqrt(1-beta2^r)/(1-beta1^r)).*Mw(3)./(sqrt(Vw(3))+epsilon);
    
    param.xi = abs(param.xi);
    param.eta = abs(param.eta);
    param.lambda = abs(param.lambda);
end


end