function [Samples] = Particle_Smoothing(QQ,data_full,data,Samples,param,target,mode,nsamples)

S = param.nParticles;

HH = param.shift;
CC = param.scale;

target.stack(:,2) = target.stack(:,2)*CC + HH;
target.stack(:,3) = target.stack(:,3)*CC;

phi_I = data_full.phi_I;
phi_C = data_full.phi_C;
phi_M = data_full.phi_M;
phi_E = data_full.phi_E;
PHI = [phi_C;phi_M;phi_E;phi_I];
PHI = log(PHI);

depth = data.depth;
depth_diff = depth(2:end) - depth(1:end-1);

d18O = data.d18O;
target.stack(:,2) = data_full.scale*target.stack(:,2) + data_full.shift;
target.stack(:,3) = data_full.scale*target.stack(:,3);
C14 = data.radiocarbon;
Age_Info = data.suggested_age;

QQ = QQ(data.index,:);

R = data_full.R;

N = length(depth);


A = cell(N,1);
W = cell(N,1);

% initialization:
[A{N},W{N}] = Proposal_PS([],[],[],d18O(N,:),C14{N},Age_Info(N,:),data_full,param,target,mode,QQ(N,:),data.ACC_MODEL,data.ACC_CONTRACTION,data.ACC_STEADY,data.ACC_EXPANSION);

% iteration:
for nn = 1:N-1
    n = N-nn;
    [A{n},W{n}] = Proposal_PS(W{n+1},A{n+1},depth_diff(n),d18O(n,:),C14{n},Age_Info(n,:),data_full,param,target,mode,QQ(n,:),data.ACC_MODEL,data.ACC_CONTRACTION,data.ACC_STEADY,data.ACC_EXPANSION);
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
    index_C = ((ZZ(n,:)==1)'.*(VV>0).*(VV<=1./1.0850)==1);
    index_M = ((ZZ(n,:)==2)'.*(VV<=1./0.9220).*(VV>1./1.0850)==1);
    index_E = ((ZZ(n,:)==3)'.*(VV>1./0.9220)==1);
    
    log_w(index_C) = log_w(index_C) + interp1(data.ACC_MODEL(:,1),data.ACC_MODEL(:,2),VV(index_C),'linear',-56) - data.ACC_CONTRACTION;
    log_w(index_M) = log_w(index_M) + interp1(data.ACC_MODEL(:,1),data.ACC_MODEL(:,2),VV(index_M),'linear',-56) - data.ACC_STEADY;
    log_w(index_E) = log_w(index_E) + interp1(data.ACC_MODEL(:,1),data.ACC_MODEL(:,2),VV(index_E),'linear',-56) - data.ACC_EXPANSION;
    
    index = index_C|index_M|index_E;
    log_w(~index) = -inf;
    
    index = (VV>1./data.lower_sedrate)|(VV<1./data.upper_sedrate);
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
index_C = ((ZZ(N-1,:)==1)'.*(VV>0).*(VV<=1./1.0850)==1);
index_M = ((ZZ(N-1,:)==2)'.*(VV<=1./0.9220).*(VV>1./1.0850)==1);
index_E = ((ZZ(N-1,:)==3)'.*(VV>1./0.9220)==1);

log_w(index_C) = log_w(index_C) + interp1(data.ACC_MODEL(:,1),data.ACC_MODEL(:,2),VV(index_C),'linear',-56) - data.ACC_CONTRACTION;
log_w(index_M) = log_w(index_M) + interp1(data.ACC_MODEL(:,1),data.ACC_MODEL(:,2),VV(index_M),'linear',-56) - data.ACC_STEADY;
log_w(index_E) = log_w(index_E) + interp1(data.ACC_MODEL(:,1),data.ACC_MODEL(:,2),VV(index_E),'linear',-56) - data.ACC_EXPANSION;

index = index_C|index_M|index_E;
log_w(~index) = -inf;

index = (VV>1./data.lower_sedrate)|(VV<1./data.upper_sedrate);
log_w(index) = -inf;

WW = exp(log_w-max(log_w,[],2));
WW = WW./sum(WW,2);

dist = cumsum(WW,2);
index = sum(rand_seed(:,n)>dist,2) + 1;

AA(N,:) = AAA(index);

% interpolation:
depth_full = data_full.depth;
Samples.ages = interp1(depth,AA,depth_full);


end