function [data,Samples,param] = getAlignment(data,data_ps,QQ,param,target,setting,MODE)

L = length(data);
data_type = setting.data_type;
nsamples_learning = setting.nSamples_learning;

Samples = struct('name',cell(L,1),'depth',cell(L,1),'ages',cell(L,1),'isoutlier',cell(L,1));
for ll = 1:L
    Samples(ll).name = data(ll).name;
    Samples(ll).depth = data(ll).depth;
end


% Learning parameters:
max_iters = 10;
threshold = 0.5;
iters = 0;
TT = ones(length(Samples),1);
while iters < max_iters && sum(TT) > 0
    iters = iters + 1;
    if strcmp(MODE,'alignment')
        disp(['#  Iteration ',num2str(iters),':']);
    end
    
    parfor ll = 1:L
        if TT(ll) == 1
            Samples(ll) = Particle_Smoothing(QQ{ll},data(ll),data_ps(ll),Samples(ll),param,target,data_type,nsamples_learning);
        end
    end
    
    parfor ll = 1:L
        if TT(ll) == 1
            [Samples(ll),QQ{ll}] = MCMC_MH(data(ll),Samples(ll),param,target,data_type,nsamples_learning);
        end
    end
    disp('   Ages are sampled.');
    
    [new_data,new_param] = updateParameters(data,Samples,param,target,setting,data_type,'align');
    
    loglik_old = getLOGLIK(data,Samples,param,target,data_type);
    loglik_new = getLOGLIK(new_data,Samples,new_param,target,data_type);
    
    new_data = UpdateParameters_R(new_data,Samples,new_param,target,TT);
    
    disp('   Parameters are updated.');
    
    if strcmp(setting.IsLearn_transition,'no')
        diff = loglik_new - loglik_old;
        tt = ['   Updating parameters makes the average log-likelihood of samples increased by [',num2str(diff'),'].'];
    elseif strcmp(setting.IsLearn_transition,'yes') == 1
        diff = mean(loglik_new) - mean(loglik_old);
        tt = ['   Updating parameters makes the average log-likelihood of samples increased by ',num2str(diff),'.'];
    end
    disp(tt);
    
    data = new_data;
    param = new_param;
    
    clear new_data;
    clear new_param;
    
    if iters > 3
        if strcmp(setting.IsLearn_transition,'no')
            TT = (abs(loglik_old-loglik_new)>threshold);
        elseif strcmp(setting.IsLearn_transition,'yes')
            if abs(diff) < threshold
                TT = zeros(length(Samples),1);
            end
        end
    end
end


% Sampling:
if strcmp(MODE,'alignment')
    nsamples = setting.nSamples;
    disp('#  Parameters are learned. Sampling algorithm is now running...');
    parfor ll = 1:L
        Samples(ll) = Particle_Smoothing(QQ{ll},data(ll),data_ps(ll),Samples(ll),param,target,data_type,nsamples);
        [Samples(ll),QQ{ll}] = MCMC_MH(data(ll),Samples(ll),param,target,data_type,nsamples);
    end
end


end