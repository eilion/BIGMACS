function [MU,NU] = learnParam_Reg(X0,X,Y,H,R,kernel,param)

M = length(X);

MU = cell(M,1);
NU = cell(M,1);

LAMBDA = param.lambda;

parfor m = 1:M
    XX0 = X0{m};
    XX = X{m};
    YY = Y{m};
    RR = LAMBDA^2*R{m};
    
    AA = XX(H{m}==0);
    YY = YY(H{m}==0);
    RR = RR(H{m}==0);
    
    N0 = size(XX0,1);
    
    K_ZZ = getCov(XX0,XX0,param,kernel) + param.lambda0^2*eye(N0);
    K_XZ = getCov(XX,XX0,param,kernel);
    K_AZ = getCov(AA,XX0,param,kernel);
    
    WW_inv = (K_ZZ+K_AZ'*(K_AZ./RR))\eye(N0);
    QQ_inv = WW_inv - K_ZZ\eye(N0);
    
    MU{m} = K_XZ*(WW_inv*(K_AZ'*(YY./RR)));
    NU{m} = param.eta^2 + param.lambda0^2 + sum((K_XZ*QQ_inv).*K_XZ,2);
end


end