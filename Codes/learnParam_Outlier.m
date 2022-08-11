function [H] = learnParam_Outlier(X,Y,MU,NU,R,param)

M = length(X);

% Detect Outliers:
H = cell(M,1);
for m = 1:M
    MEAN = MU{m};
    VAR = NU{m} + R{m};
    
    rand_seed = rand(size(X{m},1),1);
    
    DET0 = log(1-param.q) - (Y{m}-MEAN).^2./(2*sqrt(VAR).^2) - log(sqrt(VAR));
    DET1 = log(param.q/2) + log(exp(-(Y{m}-MEAN-param.d*sqrt(VAR)).^2./(2*sqrt(VAR).^2))+exp(-(Y{m}-MEAN+param.d*sqrt(VAR)).^2./(2*sqrt(VAR).^2))) - log(sqrt(VAR));
    
    AMAX = max(DET0,DET1);
    DET0 = exp(DET0-AMAX);
    DET1 = exp(DET1-AMAX);
    DET = DET1./(DET0+DET1);
    H{m} = (rand_seed<DET).*1;
    H{m}(isnan(Y{m})) = 1;
end


end