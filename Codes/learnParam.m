function [param] = learnParam(X,Y,kernel)

beta1 = 0.9;
beta2 = 0.999;
epsilon = 1e-8;
gamma = 1e-3;
max_iters = 1000;

if strcmp(kernel,'NN') == 1
    param = struct('eta',cell(1,1),'sig',cell(1,1),'lambda',cell(1,1));
    param.eta = 1;
    param.sig = [1;1];
    param.lambda = 1;
else
    param = struct('eta',cell(1,1),'xi',cell(1,1),'lambda',cell(1,1));
    param.eta = 1;
    param.xi = 1;
    param.lambda = 1;
end

if strcmp(kernel,'NN') == 1
    Mw = zeros(4,1);
    Vw = zeros(4,1);
    PDEV = zeros(4,1);
else
    Mw = zeros(3,1);
    Vw = zeros(3,1);
    PDEV = zeros(3,1);
end

M = size(X,2);

SCH = ceil(M*rand(1,max_iters));

for r = 1:max_iters
    m = SCH(r);
    XX = X{m};
    YY = Y{m};
    
    % regularlize if NN:
    if strcmp(kernel,'NN') == 1
        XX = (XX-XX(1))/(XX(end)-XX(1));
    end
    
    N = length(XX);
    rand_seed = (rand(N,1)<150/N);
    
    XX = XX(rand_seed);
    YY = YY(rand_seed);
    
    N = length(XX);
    
    KC = getCov(XX,XX,param,kernel);
    W = param.lambda^2*eye(N) + KC;
    W_inv = W\eye(N);
    alpha = W_inv*YY;
    
    if strcmp(kernel,'NN') == 1
        F = [ones(N,1),XX];
        
        SIG = diag(param.sig.^2);
        BL1 = 2*F*SIG*F';
        BL2 = 1./sqrt(1+diag(BL1));
        BL3 = BL2';
        L = BL1.*BL2.*BL3;
        
        % eta:
        WW = 2*param.eta*asin(L);
        AA = W_inv*WW*alpha;
        BB = diag(W_inv*WW*W_inv);
        PDEV(1) = sum((alpha.*AA-0.5*(1+(alpha.^2)./diag(W_inv)).*BB)./diag(W_inv));
        
        % sig:
        for t = 1:2
            BB1 = F(:,t).*F(:,t)'*2*param.sig(t);
            BB2 = diag(BB1);
            BB3 = BB2';
            
            AA = BB1.*BL2.*BL3 - 0.5*BL1.*(BL2.^3).*BB2.*BL3 - 0.5*BL1.*BL2.*(BL3.^3).*BB3;
            WW = param.eta^2./sqrt(1-L.^2).*AA;
            
            AA = W_inv*WW*alpha;
            BB = diag(W_inv*WW*W_inv);
            PDEV(1+t) = sum((alpha.*AA-0.5*(1+(alpha.^2)./diag(W_inv)).*BB)./diag(W_inv));
        end
        
        % lambda:
        WW = 2*param.lambda*eye(N);
        AA = W_inv*WW*alpha;
        BB = diag(W_inv*WW*W_inv);
        PDEV(end) = sum((alpha.*AA-0.5*(1+(alpha.^2)./diag(W_inv)).*BB)./diag(W_inv));
    else
        DIST = abs(XX-XX');
        
        if strcmp(kernel,'OU') == 1
            QQ = exp(-param.xi^2.*DIST);
        elseif strcmp(kernel,'SE') == 1
            QQ = exp(-0.5*param.xi^4.*DIST.^2);
        elseif strcmp(kernel,'M15') == 1
            QQ = (1+sqrt(3)*param.xi^2.*DIST).*exp(-sqrt(3)*param.xi^2.*DIST);
        elseif strcmp(kernel,'M25') == 1
            QQ = (1+sqrt(5)*param.xi^2.*DIST+(sqrt(5)*param.xi^2*DIST/sqrt(3)).^2).*exp(-sqrt(5)*param.xi^2.*DIST);
        end
        
        % eta:
        WW = 2*param.eta*QQ;
        AA = W_inv*WW*alpha;
        BB = diag(W_inv*WW*W_inv);
        PDEV(1) = sum((alpha.*AA-0.5*(1+(alpha.^2)./diag(W_inv)).*BB)./diag(W_inv));
        
        % xi:
        if strcmp(kernel,'OU') == 1
            WW = param.eta^2.*exp(-param.xi^2.*DIST).*(-2*param.xi.*DIST);
        elseif strcmp(kernel,'SE') == 1
            WW = param.eta^2.*exp(-0.5*param.xi^4.*DIST.^2).*(-2*param.xi^3.*DIST.^2);
        elseif strcmp(kernel,'M15') == 1
            WW = param.eta^2.*exp(-sqrt(3)*param.xi^2.*DIST).*(-2*sqrt(3)*param.xi^3.*DIST.^2);
        elseif strcmp(kernel,'M25') == 1
            WW = param.eta^2.*exp(-sqrt(5)*param.xi^2.*DIST).*(-10*param.xi^3.*DIST.^2)/3.*(1+sqrt(5)*param.xi^2.*DIST);
        end
        AA = W_inv*WW*alpha;
        BB = diag(W_inv*WW*W_inv);
        PDEV(2) = sum((alpha.*AA-0.5*(1+(alpha.^2)./diag(W_inv)).*BB)./diag(W_inv));
        
        % lambda:
        WW = 2*param.lambda*eye(N);
        AA = W_inv*WW*alpha;
        BB = diag(W_inv*WW*W_inv);
        PDEV(end) = sum((alpha.*AA-0.5*(1+(alpha.^2)./diag(W_inv)).*BB)./diag(W_inv));
    end
    
    Mw = beta1*Mw - (1-beta1)*PDEV;
    Vw = beta2*Vw + (1-beta2)*PDEV.*PDEV;
    
    if strcmp(kernel,'NN') == 1
        param.eta = param.eta - (gamma*sqrt(1-beta2^r)/(1-beta1^r)).*Mw(1)./(sqrt(Vw(1))+epsilon);
        param.sig = param.sig - (gamma*sqrt(1-beta2^r)/(1-beta1^r)).*Mw(2:end-1)./(sqrt(Vw(2:end-1))+epsilon);
        param.lambda = param.lambda - (gamma*sqrt(1-beta2^r)/(1-beta1^r)).*Mw(end)./(sqrt(Vw(end))+epsilon);
        param.sig = abs(param.sig);
    else
        param.eta = param.eta - (gamma*sqrt(1-beta2^r)/(1-beta1^r)).*Mw(1)./(sqrt(Vw(1))+epsilon);
        param.xi = param.xi - (gamma*sqrt(1-beta2^r)/(1-beta1^r)).*Mw(2)./(sqrt(Vw(2))+epsilon);
        param.lambda = param.lambda - (gamma*sqrt(1-beta2^r)/(1-beta1^r)).*Mw(3)./(sqrt(Vw(3))+epsilon);
        param.xi = abs(param.xi);
    end
    param.eta = abs(param.eta);
    param.lambda = abs(param.lambda);
end


end