function [new_stack,rcd_stack] = Stack_Updater_2(Samples,data,initial_stack,stack,core_param,kernel)

L = length(Samples);
TT = size(stack,1);

MEAN = stack(:,2);
MEAN_V = mean(MEAN);

param = struct('eta',cell(L,1),'xi',cell(L,1),'lambda',cell(L,1));
rcd_stack = struct('name',cell(L,1),'age',cell(L,1),'mean',cell(L,1),'stdv',cell(L,1));

AAA = NaN*ones(TT,100*L);
BBB = NaN*ones(TT,100*L);

for ll = 1:L
    Age = Samples(ll).ages;
    HH = Samples(ll).isoutlier;
    d18O = (data(ll).d18O-core_param(ll).shift)/core_param(ll).scale;
    
    M = size(Age,2);
    
    X = cell(M,1);
    Y = cell(M,1);
    
    for m = 1:M
        A = Age(:,m);
        V = d18O;
        H = HH{m};
        A = repmat(A,[1,size(V,2)]);
        
        index = (isnan(V)==0);
        
        A = A(index);
        V = V(index);
        H = H(index);
        
        X{m} = A(H==0);
        Y{m} = V(H==0);
    end
    
    param(ll) = learnParam(X,Y,kernel);
    
    % construct a record-specific stack:
    rcd_stack(ll).name = data(ll).name;
    
    AA = NaN*ones(TT,M);
    BB = NaN*ones(TT,M);
    ST = zeros(1,M);
    ED = zeros(1,M);
    
    for m = 1:M
        
        AA(:,m) = initial_stack(:,2);
        BB(:,m) = initial_stack(:,3).^2;
        
        [~,order] = sort(abs(stack(:,1)-min(X{m})),'ascend');
        st = order(1);
        
        [~,order] = sort(abs(stack(:,1)-max(X{m})),'ascend');
        ed = order(1);
        
        ST(m) = st;
        ED(m) = ed;
        
        XX = stack(st:ed,1);
        
        T = X{m};
        V = Y{m} - MEAN_V;
        
        KV = param(ll).lambda^2*eye(length(T)) + getCov(T,T,param(ll),kernel);
        KXV = getCov(XX,T,param(ll),kernel);
        KXX = getCov(XX,XX,param(ll),kernel);
        
        Mu = KXV*(KV\V) + MEAN_V;
        Sigma = param(ll).lambda^2*eye(length(XX)) + KXX - KXV*(KV\KXV');
        
        AA(st:ed,m) = Mu;
        BB(st:ed,m) = diag(Sigma);
    end
    
    DD = mean(AA,2);
    EE = mean(BB+(AA-DD).^2,2);
    
    st = max(ST);
    ed = min(ED);
    
    rcd_stack(ll).age = stack(st:ed,1);
    rcd_stack(ll).mean = DD(st:ed);
    rcd_stack(ll).stdv = sqrt(EE(st:ed));
    
    rand_seed = ceil(M*rand(1,100));
    
    for k = 1:100
        AAA(ST(rand_seed(k)):ED(rand_seed(k)),(ll-1)*100+k) = DD(ST(rand_seed(k)):ED(rand_seed(k)));
        BBB(ST(rand_seed(k)):ED(rand_seed(k)),(ll-1)*100+k) = EE(ST(rand_seed(k)):ED(rand_seed(k)));
    end
end

% construct the local stack:
A = zeros(TT,1);
B = zeros(TT,1);
C = zeros(TT,1);

for ll = 1:100*L
    index = (isnan(AAA(:,ll))==0);
    A(index) = A(index) + sum(AAA(index,ll)./BBB(index,ll),2);
    B(index) = B(index) + sum(1./BBB(index,ll),2);
    C(index) = C(index) + 1;
end

index = (C>0);

new_stack = stack;
new_stack(index,2) = A(index)./B(index);
new_stack(index,3) = sqrt(C(index)./B(index));


end

