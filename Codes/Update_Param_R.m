function [core_param] = Update_Param_R(data,Samples,core_param,stack,TTT)

L = length(data);

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
                
                param = learnParam(X,Y,'NN');
                
                %{
            % K = min(10,M);
            K = 1;
            IND = zeros(max_iters,K);
            for k = 1:max_iters
                IND(k,:) = randperm(M,K);
            end
                %}
                MU = zeros(size(stack,1),1);
                
                for m = 1:M
                    % XX = cat(1,X{IND(m,:)});
                    % RR = cat(1,Y{IND(m,:)});
                    XX = cat(1,X{m});
                    RR = cat(1,Y{m});
                    
                    rand_seed = (rand(length(XX),1)<150/length(XX));
                    XX = XX(rand_seed);
                    RR = RR(rand_seed);
                    
                    [~,order] = sort(XX,'ascend');
                    XX = XX(order);
                    RR = RR(order);
                    
                    T = core_param(ll).R(:,1);
                    
                    TT = (T-XX(1))/(XX(end)-XX(1));
                    XX = (XX-XX(1))/(XX(end)-XX(1));
                    
                    COV_TX = getCov(TT,XX,param,'NN');
                    COV_XX = getCov(XX,XX,param,'NN') + param.lambda^2*eye(size(XX,1));
                    
                    MU = MU + COV_TX*(COV_XX\RR);
                end
                
                core_param(ll).R(:,2) = exp(MU/M);
            end
        end
    end
end


end