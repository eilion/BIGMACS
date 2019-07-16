function [core_param] = Update_Param_R_2(data,Samples,global_param,core_param,stack,TTT)

L = length(data);

h = global_param.h;

for ll = 1:L
    [~,~,~,~,~,islearn_average_sed_rate,~,~] = getSetting_Core(data(ll).name);
    if TTT(ll) == 1
        if strcmp(islearn_average_sed_rate,'no') == 0
            if strcmp(islearn_average_sed_rate,'constant') == 1
                X = Samples(ll).ages;
                D = data(ll).depth;
                AA = (X(end,:)-X(1,:))./(D(end)-D(1));
                core_param(ll).R(:,2) = mean(AA)*ones(size(stack,1),1);
            elseif strcmp(islearn_average_sed_rate,'adaptive') == 1
                M = size(Samples(ll).ages,2);
                
                Age = Samples(ll).ages;
                depth = data(ll).depth;
                
                index = (isnan(Age(:,1))==0);
                A = Age(index,:);
                depth = depth(index);
                
                R = (A(2:end,:)-A(1:end-1,:))./(depth(2:end)-depth(1:end-1));
                A = A(2:end,:);
                
                X = cell(M,1);
                Y = cell(M,1);
                
                for m = 1:M
                    index = ((isnan(R(:,m))==0)&(isinf(R(:,m))==0));
                    X{m} = A(index,m);
                    Y{m} = log(R(index,m));
                end
                
                T = core_param(ll).R(:,1);
                MU = zeros(size(T,1),1);
                
                for m = 1:M
                    AA = - 0.5*((T-X{m}')/h).^2 - log(h) - 0.5*log(2*pi);
                    amax = max(AA,[],2);
                    BB = amax + log(sum(exp(AA-amax),2));
                    
                    MU = MU + sum(exp(AA-BB).*Y{m}',2);
                end
                
                core_param(ll).R(:,2) = exp(MU/M);
            end
        end
    end
end


end