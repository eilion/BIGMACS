function [QQ] = getProposal(data_full,data,param,target,data_type)

S = 500;
nsamples = 100;

HH = param.shift;
CC = param.scale;

phi_I = data.phi_I;
phi_C = data.phi_C;
phi_M = data.phi_M;
phi_E = data.phi_E;
PHI = [phi_C;phi_M;phi_E;phi_I];
PHI = log(PHI);

depth = data.depth;
depth_diff = depth(2:end) - depth(1:end-1);
d18O = data.d18O;

target.stack(:,2) = target.stack(:,2)*CC + HH;
target.stack(:,3) = target.stack(:,3)*CC;

target.stack(:,2) = data.scale*target.stack(:,2) + data.shift;
target.stack(:,3) = data.scale*target.stack(:,3);
C14 = data.radiocarbon;
Age_Info = data.suggested_age;
Age_Info = [Age_Info,zeros(size(Age_Info,1),2)];

for n = 1:size(Age_Info,1)
    if ~isnan(Age_Info(n,1))
        if Age_Info(n,3) == 0
            Age_Info(n,4) = max(data.min,Age_Info(n,1)-3*Age_Info(n,2));
        else
            Age_Info(n,4) = max(data.min,Age_Info(n,1)-Age_Info(n,2));
        end
    else
        if n == 1
            Age_Info(n,4) = data.min;
        else
            Age_Info(n,4) = Age_Info(n-1,4);
        end
    end
end

for nn = 1:size(Age_Info,1)
    n = size(Age_Info,1) - nn + 1;
    if ~isnan(Age_Info(n,1))
        if Age_Info(n,3) == 0
            Age_Info(n,5) = min(data.max,Age_Info(n,1)+3*Age_Info(n,2));
        else
            Age_Info(n,5) = min(data.max,Age_Info(n,1)+Age_Info(n,2));
        end
    else
        if n == size(Age_Info,1)
            Age_Info(n,5) = data.max;
        else
            Age_Info(n,5) = Age_Info(n+1,5);
        end
    end
end

R = data.R;

N = length(depth);


A = cell(N,1);
W = cell(N,1);

% initialization:
[A{N},W{N}] = Proposal_init([],[],[],d18O(N,:),C14{N},Age_Info(N,:),data,param,S,target,data_type,N);

% iteration:
for nn = 1:N-1
    n = N-nn;
    [A{n},W{n}] = Proposal_init(W{n+1},A{n+1},depth_diff(n),d18O(n,:),C14{n},Age_Info(n,:),data,param,S,target,data_type,n);
end


% backward sampling:
ZZZ = repelem([1,2,3],S);
AA = zeros(N,nsamples);
ZZ = zeros(N,nsamples);

AAA = reshape(A{1}',[1,3*S]);
WW = reshape(W{1}',[1,3*S]);
WW = exp(WW-max(WW));
WW = WW./sum(WW);

index = randsample(3*S,nsamples,true,WW);
AA(1,:) = AAA(index);
ZZ(1,:) = ZZZ(index);

rand_seed = rand(nsamples,N-1);
for n = 1:N-2
    AAA = reshape(A{n+1}',[1,3*S]);
    
    log_w = zeros(nsamples,3*S);
    log_w(:,1:S) = W{n+1}(1,:) + PHI(1,ZZ(n,:))';
    log_w(:,S+1:2*S) = W{n+1}(2,:) + PHI(2,ZZ(n,:))';
    log_w(:,2*S+1:3*S) = W{n+1}(3,:) + PHI(3,ZZ(n,:))';
    
    RR = interp1(R(:,1),R(:,2),AAA);
    RR = repmat(RR,[nsamples,1]);
    VV = (AAA-AA(n,:)')./(RR.*depth_diff(n));
    index_C = ((ZZ(n,:)==1)'.*(VV>0).*(VV<0.9220)==1);
    index_M = ((ZZ(n,:)==2)'.*(VV>=0.9220).*(VV<1.0850)==1);
    index_E = ((ZZ(n,:)==3)'.*(VV>=1.0850)==1);
    
    log_w(index_C) = log_w(index_C) + interp1(data.ACC_MODEL(:,1),data.ACC_MODEL(:,2),VV(index_C),'linear',-56) - data.ACC_CONTRACTION;
    log_w(index_M) = log_w(index_M) + interp1(data.ACC_MODEL(:,1),data.ACC_MODEL(:,2),VV(index_M),'linear',-56) - data.ACC_STEADY;
    log_w(index_E) = log_w(index_E) + interp1(data.ACC_MODEL(:,1),data.ACC_MODEL(:,2),VV(index_E),'linear',-56) - data.ACC_EXPANSION;
    
    index = index_C|index_M|index_E;
    log_w(~index) = -inf;
    
    index = (VV<data.lower_sedrate)|(VV>data.upper_sedrate);
    log_w(index) = -inf;
    
    WW = exp(log_w-max(log_w,[],2));
    WW = WW./sum(WW,2);
    
    dist = cumsum(WW,2);
    index = sum(rand_seed(:,n)>dist,2) + 1;
    
    AA(n+1,:) = AAA(index);
    ZZ(n+1,:) = ZZZ(index);
end

AAA = A{N};

log_w = W{N} + PHI(4,ZZ(n,:))';

RR = interp1(R(:,1),R(:,2),AAA);
RR = repmat(RR,[nsamples,1]);
VV = (AAA-AA(N-1,:)')./(RR.*depth_diff(N-1));
index_C = ((ZZ(N-1,:)==1)'.*(VV>0).*(VV<0.9220)==1);
index_M = ((ZZ(N-1,:)==2)'.*(VV>=0.9220).*(VV<1.0850)==1);
index_E = ((ZZ(N-1,:)==3)'.*(VV>=1.0850)==1);
log_w(index_C) = log_w(index_C) + interp1(data.ACC_MODEL(:,1),data.ACC_MODEL(:,2),VV(index_C),'linear',-56) - data.ACC_CONTRACTION;
log_w(index_M) = log_w(index_M) + interp1(data.ACC_MODEL(:,1),data.ACC_MODEL(:,2),VV(index_M),'linear',-56) - data.ACC_STEADY;
log_w(index_E) = log_w(index_E) + interp1(data.ACC_MODEL(:,1),data.ACC_MODEL(:,2),VV(index_E),'linear',-56) - data.ACC_EXPANSION;

index = index_C|index_M|index_E;
log_w(~index) = -inf;

index = (VV<data.lower_sedrate)|(VV>data.upper_sedrate);
log_w(index) = -inf;

WW = exp(log_w-max(log_w,[],2));
WW = WW./sum(WW,2);

dist = cumsum(WW,2);
index = sum(rand_seed(:,n)>dist,2) + 1;

AA(N,:) = AAA(index);

% interpolation:
depth_full = data_full.depth;
N = length(depth_full);
AAA = interp1(depth,AA,depth_full);

QQ = zeros(N,2);
QQ(:,1) = mean(AAA,2);
QQ(:,2) = sqrt(var(AAA,1,2)) + 1.285*data.PTCL_BW;


end