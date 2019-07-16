function [K] = getCov(X1,X2,param,kernel)

if strcmp(kernel,'NN') == 1
    N1 = size(X1,1);
    N2 = size(X2,1);
    XX = [ones(N1,1),X1];
    YY = [ones(N2,1),X2];
    
    SIG = diag(param.sig.^2);
    AA = XX*SIG*YY';
    BB = diag(XX*SIG*XX');
    CC = diag(YY*SIG*YY');
    
    K = param.eta^2*asin(2*AA./sqrt((1+2*BB).*(1+2*CC')));
else
    DIST = abs(X1-X2');
    if strcmp(kernel,'OU') == 1
        K = param.eta^2*exp(-param.xi^2*DIST);
    elseif strcmp(kernel,'SE') == 1
        K = param.eta^2*exp(-0.5*param.xi^4*DIST.^2);
    elseif strcmp(kernel,'M15') == 1
        K = param.eta^2*(1+sqrt(3)*param.xi^2*DIST).*exp(-sqrt(3)*param.xi^2*DIST);
    elseif strcmp(kernel,'M25') == 1
        K = param.eta^2*(1+sqrt(5)*param.xi^2*DIST+(sqrt(5)*param.xi^2*DIST/sqrt(3)).^2).*exp(-sqrt(5)*param.xi^2*DIST);
    end
end


end