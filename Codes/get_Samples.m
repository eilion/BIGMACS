function [Samples,QQ] = get_Samples(Samples,QQ,data,core_param,global_param,cal_curve,stack,mode,nsamples,TT)

L = length(data);

parfor ll = 1:L
    if TT(ll) == 1
        [Samples(ll),~] = PS(QQ{ll},data(ll),Samples(ll),core_param(ll),global_param,stack,cal_curve,mode,nsamples);
        [Samples(ll),QQ{ll}] = MCMC_MH(data(ll),Samples(ll),core_param(ll),global_param,stack,cal_curve,mode,nsamples);
    end
end


end