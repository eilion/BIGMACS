function [Samples,QQ] = PS(QQ,data,Samples,core_param,global_param,stack,cal_curve,mode,nsamples)

S = global_param.nParticles;

alpha = global_param.alpha;
beta = global_param.beta;

phi_I = core_param.phi_I;
phi_C = core_param.phi_C;
phi_M = core_param.phi_M;
phi_E = core_param.phi_E;
PHI = [phi_C;phi_M;phi_E;phi_I];
PHI = log(PHI);

depth = data.depth;
depth_diff = depth(2:end) - depth(1:end-1);
% d18O = (data.d18O-core_param.shift)/core_param.scale;
d18O = data.d18O;
stack(:,2) = core_param.scale*stack(:,2) + core_param.shift;
stack(:,3) = core_param.scale*stack(:,3);
C14 = data.radiocarbon;
Age_Info = data.suggested_age;

R = core_param.R;

N = length(depth);


A = cell(N,1);
W = cell(N,1);

% initialization:
[A{N},W{N}] = Proposal_5([],[],[],d18O(N,:),C14{N},Age_Info(N,:),core_param,global_param,stack,cal_curve,mode,QQ(N,:));

% iteration:
for nn = 1:N-1
    n = N-nn;
    [A{n},W{n}] = Proposal_5(W{n+1},A{n+1},depth_diff(n),d18O(n,:),C14{n},Age_Info(n,:),core_param,global_param,stack,cal_curve,mode,QQ(n,:));
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
    log_w(index_C) = log_w(index_C) + (alpha-1)*log(VV(index_C)) - beta*VV(index_C) - log(RR(index_C)) - log(gamcdf(0.9220,alpha,1/beta));
    log_w(index_M) = log_w(index_M) + (alpha-1)*log(VV(index_M)) - beta*VV(index_M) - log(RR(index_M)) - log(gamcdf(1.0850,alpha,1/beta)-gamcdf(0.9220,alpha,1/beta));
    log_w(index_E) = log_w(index_E) + (alpha-1)*log(VV(index_E)) - beta*VV(index_E) - log(RR(index_E)) - log(1-gamcdf(1.0850,alpha,1/beta));
    
    index = index_C|index_M|index_E;
    log_w(~index) = -inf;
    
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
log_w(index_C) = log_w(index_C) + (alpha-1)*log(VV(index_C)) - beta*VV(index_C) - log(RR(index_C)) - log(gamcdf(0.9220,alpha,1/beta));
log_w(index_M) = log_w(index_M) + (alpha-1)*log(VV(index_M)) - beta*VV(index_M) - log(RR(index_M)) - log(gamcdf(1.0850,alpha,1/beta)-gamcdf(0.9220,alpha,1/beta));
log_w(index_E) = log_w(index_E) + (alpha-1)*log(VV(index_E)) - beta*VV(index_E) - log(RR(index_E)) - log(1-gamcdf(1.0850,alpha,1/beta));

index = index_C|index_M|index_E;
log_w(~index) = -inf;

WW = exp(log_w-max(log_w,[],2));
WW = WW./sum(WW,2);

dist = cumsum(WW,2);
index = sum(rand_seed(:,n)>dist,2) + 1;

AA(N,:) = AAA(index);

Samples.ages = AA;

%{
QQ = zeros(N,2);
QQ(:,1) = mean(AA,2);
QQ(:,2) = var(AA,1,2) + 1;
%}
QQ = zeros(N,S);
index = ceil(nsamples*rand(N,S));
for n = 1:N
    QQ(n,:) = AA(n,index(n,:));
end


end