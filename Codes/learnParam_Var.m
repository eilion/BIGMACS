function [param,R] = learnParam_Var(X,Y,H,MU,NU,param,variance,r)

M = length(X);
R = cell(M,1);

if strcmp(variance,'homoscedastic') == 1 || r == 1
    for m = 1:M
        N = size(X{m},1);
        R{m} = ones(N,1);
    end
elseif strcmp(variance,'heteroscedastic') == 1
    % Decide K:
    KK_CAND = zeros(M,1);
    
    parfor m = 1:M
        XX = X{m}(H{m}==0);
        YY = Y{m}(H{m}==0);
        
        MEAN = MU{m}(H{m}==0);
        VAR = NU{m}(H{m}==0);
        
        NN = size(XX,1);
        
        QQ = abs(XX-XX');
        PP = sort(QQ,2,'ascend');
        
        KK = (ceil(500/NN):1:max(30,ceil(300/NN)));
        
        LOGLIK = zeros(length(KK),1);
        
        for k = 1:length(KK)
            BW = PP(:,ceil(NN*KK(k)/100));
            
            WW = - 0.5*QQ.^2./(BW'.^2) - log(BW');
            WW(abs(QQ)<1e-12) = -inf;
            WW = WW - max(WW,[],2);
            WW = exp(WW);
            LL = WW./sum(WW,2);
            
            LA = sum(((YY'-MEAN').^2+VAR').*LL,2);
            
            LOGLIK(k) = - 0.5*sum((YY-MEAN).^2.*(VAR+LA).^(-1)) - 0.5*sum(log(VAR+LA));
        end
        
        [~,amax] = max(LOGLIK);
        KK_CAND(m) = KK(amax);
    end
    
    param.K = KK_CAND;
    
    
    clear LA;
    
    KK = param.K;
    
    % Update R:
    parfor m = 1:M
        XX = X{m}(H{m}==0);
        YY = Y{m}(H{m}==0);
        
        MEAN = MU{m}(H{m}==0);
        VAR = NU{m}(H{m}==0);
        
        NN = size(XX,1);
        
        QQ = abs(XX-XX');
        PP = sort(QQ,2,'ascend');
        BW = PP(:,ceil(NN*KK(m)/100));
        
        QQ = abs(X{m}-XX');
        
        WW = - 0.5*QQ.^2./(BW'.^2) - log(BW');
        WW = WW - max(WW,[],2);
        WW = exp(WW);
        LL = WW./sum(WW,2);
        
        R{m} = sum(((YY'-MEAN').^2+VAR').*LL,2);
    end
    
    param.lambda = 1;
end


end