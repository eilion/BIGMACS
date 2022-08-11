function [data] = UpdateParameters_R(data,Samples,param,target,TTT)

L = length(data);

h = param.h;

for ll = 1:L
    if TTT(ll) == 1
        if strcmp(data(ll).islearn_average_sed_rate,'no') == 0
            if strcmp(data(ll).islearn_average_sed_rate,'constant') == 1
                X = Samples(ll).ages;
                D = data(ll).depth;
                AA = (X(end,:)-X(1,:))./(D(end)-D(1));
                data(ll).R(:,2) = mean(AA)*ones(size(target.stack,1),1);
            elseif strcmp(data(ll).islearn_average_sed_rate,'adaptive') == 1
                M = size(Samples(ll).ages,2);
                
                Age = Samples(ll).ages;
                depth = data(ll).depth;
                
                index = ~isnan(Age(:,1));
                A = Age(index,:);
                depth = depth(index);
                
                R = (A(2:end,:)-A(1:end-1,:))./(depth(2:end)-depth(1:end-1));
                A = A(2:end,:);
                
                X = cell(M,1);
                Y = cell(M,1);
                
                for m = 1:M
                    index = (~isnan(R(:,m))&~isinf(R(:,m)));
                    X{m} = A(index,m);
                    Y{m} = log(R(index,m));
                end
                
                T = linspace(data(ll).R(1,1),data(ll).R(end,1),1000)';
                % T = data(ll).R(:,1);
                MU = zeros(size(T,1),1);
                
                for m = 1:M
                    AA = - 0.5*((T-X{m}')/h).^2 - log(h) - 0.5*log(2*pi);
                    amax = max(AA,[],2);
                    BB = amax + log(sum(exp(AA-amax),2));
                    
                    MU = MU + sum(exp(AA-BB).*Y{m}',2);
                end
                
                data(ll).R(:,2) = interp1(T,exp(MU/M),data(ll).R(:,1));
            end
        end
    end
end


end