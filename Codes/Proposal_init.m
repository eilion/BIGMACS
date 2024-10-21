function [SAM_A,WW] = Proposal_init(WWB,AB,depth_diff,d18O,C14_Table,Age_Info,data,param,S,target,mode,nn)

stack = target.stack;
cal_curve = target.cal_curve;

a_d18O = param.a_d18O;
b_d18O = param.b_d18O;
log_tail = gammaln(a_d18O+0.5) - gammaln(a_d18O) - 0.5*log(2*pi*b_d18O);

a = param.a_C14;
b = param.b_C14;

R = data.R;

phi_I = data.phi_I;
phi_C = data.phi_C;
phi_M = data.phi_M;
phi_E = data.phi_E;
PHI = [phi_C;phi_M;phi_E;phi_I];
PHI = log(PHI);

% data.lower_sedrate = 0;
% data.upper_sedrate = inf;

if isempty(WWB) == 1
    age_ed = Age_Info(5);
    age_st = Age_Info(4);
    
    if ~isempty(data.initial_age)
        ID = (abs(data.initial_age(:,1)-data.depth(nn))<1e-6);
        if sum(ID) > 0
            ZZ = data.initial_age(ID,2);
            if ZZ(1) > age_st && ZZ(1) < age_ed
                SAM_A = ZZ(1)*ones(1,S);
            else
                SAM_A = (age_ed-age_st).*rand(1,S) + age_st;
            end
        else
            age_range = age_ed - age_st;
            % Sample ages:
            SAM_A = age_range*rand(1,S) + age_st;
        end
    else
        age_range = age_ed - age_st;
        % Sample ages:
        SAM_A = age_range*rand(1,S) + age_st;
    end
    
    % Compute weights:
    MargLik = zeros(1,S);
    if strcmp(mode,'C14') == 0 && ~isnan(d18O(1))
        mu = interp1(stack(:,1),stack(:,2),SAM_A);
        sig = interp1(stack(:,1),stack(:,3),SAM_A);
        N = sum(~isnan(d18O));
        for n = 1:N
            % AA = (1-q)*exp(-(d18O(n)-mu).^2./(2*sig.^2))./sqrt(2*pi*sig.^2) + q/2*exp(-(d18O(n)-mu-d*sig).^2./(2*sig.^2))./sqrt(2*pi*sig.^2) + q/2*exp(-(d18O(n)-mu+d*sig).^2./(2*sig.^2))./sqrt(2*pi*sig.^2);
            % MargLik = MargLik + log(AA);
            MargLik = MargLik - (a_d18O+0.5)*log(1+(d18O(n)-mu).^2./(2*b_d18O*sig.^2)) - log(sig) + log_tail;
        end
    end
    if strcmp(mode,'d18O') == 0 && ~isempty(C14_Table)
        MargLik = MargLik + EP_Y(SAM_A,C14_Table,cal_curve,a,b);
    end
    
    if ~isnan(Age_Info(1))
        if Age_Info(3) == 0
            MargLik = MargLik - (SAM_A-Age_Info(1)).^2./(2*Age_Info(2).^2);
        elseif Age_Info(3) == 1
            index = (SAM_A>=Age_Info(1)-Age_Info(2))&(SAM_A<=Age_Info(1)+Age_Info(2));
            MargLik(~index) = -inf;
        end
    end
    
    MargLik(SAM_A<age_st|SAM_A>age_ed) = -inf;
    MargLik(isnan(MargLik)) = -inf;
    
    if ~isnan(data.max)
        MargLik(SAM_A>data.max) = -inf;
    end
    if ~isnan(data.min)
        MargLik(SAM_A<data.min) = -inf;
    end
    
    WW = MargLik;
else
    RR = interp1(R(:,1),R(:,2),min(AB(~isinf(WWB))));
    age_st = max(Age_Info(4),min(AB(~isinf(WWB)))-1./data.lower_sedrate*RR*depth_diff);
    RR = interp1(R(:,1),R(:,2),max(AB(~isinf(WWB))));
    age_ed = max(age_st,min(Age_Info(5),max(AB(~isinf(WWB)))-1./data.upper_sedrate*RR*depth_diff));
    
    if size(WWB,1) == 1
        index = ~isinf(WWB);
        AB_TABLE = AB(index);
        WB_TABLE = WWB(index);
        
        % AB_TABLE = AB(1:min(200,S));
        % WB_TABLE = WWB(1:min(200,S));
        % AB_TABLE = AB;
        % WB_TABLE = WWB;
        
        if ~isempty(data.initial_age)
            ID = (abs(data.initial_age(:,1)-data.depth(nn))<1e-6);
            if sum(ID) > 0
                ZZ = data.initial_age(ID,2);
                if ZZ(1) > age_st && ZZ(1) < age_ed
                    SAM_A = ZZ(1)*ones(3,S);
                else
                    SAM_A = (age_ed-age_st).*rand(3,S) + age_st;
                end
            else
                % Sample ages:
                SAM_A = (age_ed-age_st).*rand(3,S) + age_st;
            end
        else
            % Sample ages:
            SAM_A = (age_ed-age_st).*rand(3,S) + age_st;
        end
        
        % Compute weights:
        WW = zeros(3,S);
        
        RR = interp1(R(:,1),R(:,2),AB_TABLE');
        RR = repmat(RR,[1,S]);
        
        % Z == 1:
        MargLik = PHI(4,1) + WB_TABLE';
        
        VV = (AB_TABLE'-SAM_A(1,:))./(RR.*depth_diff);
        
        MargLik = MargLik + interp1(data.ACC_MODEL(:,1),data.ACC_MODEL(:,2),VV,'linear',-56) - data.ACC_CONTRACTION;
        index = (VV<=0)|(VV>1./1.0850);
        MargLik(index) = -inf;
        
        index = (VV>1./data.lower_sedrate+0.001)|(VV<1./data.upper_sedrate-0.001);
        MargLik(index) = -inf;
        
        AMAX = max(MargLik);
        MargLik = AMAX + log(sum(exp(MargLik-AMAX)));
        
        if strcmp(mode,'C14') == 0 && ~isnan(d18O(1))
            % MargLik = MargLik + log((1-q)*exp(EP_V(SAM_A(1,:),zeros(1,S),d18O,stack))+q*exp(EP_V(SAM_A(1,:),ones(1,S),d18O,stack)));
            mu = interp1(stack(:,1),stack(:,2),SAM_A(1,:));
            sig = interp1(stack(:,1),stack(:,3),SAM_A(1,:));
            N = sum(~isnan(d18O));
            for n = 1:N
                % AA = (1-q)*exp(-(d18O(n)-mu).^2./(2*sig.^2))./sqrt(2*pi*sig.^2) + q/2*exp(-(d18O(n)-mu-d*sig).^2./(2*sig.^2))./sqrt(2*pi*sig.^2) + q/2*exp(-(d18O(n)-mu+d*sig).^2./(2*sig.^2))./sqrt(2*pi*sig.^2);
                % MargLik = MargLik + log(AA);
                MargLik = MargLik - (a_d18O+0.5)*log(1+(d18O(n)-mu).^2./(2*b_d18O*sig.^2)) - log(sig) + log_tail;
            end
        end
        if strcmp(mode,'d18O') == 0 && ~isempty(C14_Table)
            MargLik = MargLik + EP_Y(SAM_A(1,:),C14_Table,cal_curve,a,b);
        end
        
        if ~isnan(Age_Info(1))
            if Age_Info(3) == 0
                MargLik = MargLik - (SAM_A(1,:)-Age_Info(1)).^2./(2*Age_Info(2).^2);
            elseif Age_Info(3) == 1
                index = (SAM_A(1,:)>=Age_Info(1)-Age_Info(2))&(SAM_A(1,:)<=Age_Info(1)+Age_Info(2));
                MargLik(~index) = -inf;
            end
        end
        
        MargLik(SAM_A(1,:)<age_st) = -inf;
        MargLik(isnan(MargLik)) = -inf;
        
        if ~isnan(data.max)
            MargLik(SAM_A(1,:)>data.max) = -inf;
        end
        if ~isnan(data.min)
            MargLik(SAM_A(1,:)<data.min) = -inf;
        end
        
        WW(1,:) = MargLik;
        
        % Z == 2:
        MargLik = PHI(4,2) + WB_TABLE';
        
        VV = (AB_TABLE'-SAM_A(2,:))./(RR.*depth_diff);
        
        MargLik = MargLik + interp1(data.ACC_MODEL(:,1),data.ACC_MODEL(:,2),VV,'linear',-56) - data.ACC_STEADY;
        index = (VV<=1./1.0850)|(VV>1./0.9220);
        MargLik(index) = -inf;
        
        index = (VV>1./data.lower_sedrate+0.001)|(VV<1./data.upper_sedrate-0.001);
        MargLik(index) = -inf;
        
        AMAX = max(MargLik);
        MargLik = AMAX + log(sum(exp(MargLik-AMAX)));
        
        if strcmp(mode,'C14') == 0 && ~isnan(d18O(1))
            % MargLik = MargLik + log((1-q)*exp(EP_V(SAM_A(2,:),zeros(1,S),d18O,stack))+q*exp(EP_V(SAM_A(2,:),ones(1,S),d18O,stack)));
            mu = interp1(stack(:,1),stack(:,2),SAM_A(2,:));
            sig = interp1(stack(:,1),stack(:,3),SAM_A(2,:));
            N = sum(~isnan(d18O));
            for n = 1:N
                % AA = (1-q)*exp(-(d18O(n)-mu).^2./(2*sig.^2))./sqrt(2*pi*sig.^2) + q/2*exp(-(d18O(n)-mu-d*sig).^2./(2*sig.^2))./sqrt(2*pi*sig.^2) + q/2*exp(-(d18O(n)-mu+d*sig).^2./(2*sig.^2))./sqrt(2*pi*sig.^2);
                % MargLik = MargLik + log(AA);
                MargLik = MargLik - (a_d18O+0.5)*log(1+(d18O(n)-mu).^2./(2*b_d18O*sig.^2)) - log(sig) + log_tail;
            end
        end
        if strcmp(mode,'d18O') == 0 && ~isempty(C14_Table)
            MargLik = MargLik + EP_Y(SAM_A(2,:),C14_Table,cal_curve,a,b);
        end
        
        if ~isnan(Age_Info(1))
            if Age_Info(3) == 0
                MargLik = MargLik - (SAM_A(2,:)-Age_Info(1)).^2./(2*Age_Info(2).^2);
            elseif Age_Info(3) == 1
                index = (SAM_A(2,:)>=Age_Info(1)-Age_Info(2))&(SAM_A(2,:)<=Age_Info(1)+Age_Info(2));
                MargLik(~index) = -inf;
            end
        end
        
        MargLik(SAM_A(2,:)<age_st) = -inf;
        MargLik(isnan(MargLik)) = -inf;
        
        if ~isnan(data.max)
            MargLik(SAM_A(2,:)>data.max) = -inf;
        end
        if ~isnan(data.min)
            MargLik(SAM_A(2,:)<data.min) = -inf;
        end
        
        WW(2,:) = MargLik;
        
        % Z == 3:
        MargLik = PHI(4,3) + WB_TABLE';
        
        VV = (AB_TABLE'-SAM_A(3,:))./(RR.*depth_diff);
        
        MargLik = MargLik + interp1(data.ACC_MODEL(:,1),data.ACC_MODEL(:,2),VV,'linear',-56) - data.ACC_EXPANSION;
        index = (VV<=1./0.9220);
        MargLik(index) = -inf;
        
        index = (VV>1./data.lower_sedrate+0.001)|(VV<1./data.upper_sedrate-0.001);
        MargLik(index) = -inf;
        
        AMAX = max(MargLik);
        MargLik = AMAX + log(sum(exp(MargLik-AMAX)));
        
        if strcmp(mode,'C14') == 0 && ~isnan(d18O(1))
            % MargLik = MargLik + log((1-q)*exp(EP_V(SAM_A(3,:),zeros(1,S),d18O,stack))+q*exp(EP_V(SAM_A(3,:),ones(1,S),d18O,stack)));
            mu = interp1(stack(:,1),stack(:,2),SAM_A(3,:));
            sig = interp1(stack(:,1),stack(:,3),SAM_A(3,:));
            N = sum(~isnan(d18O));
            for n = 1:N
                % AA = (1-q)*exp(-(d18O(n)-mu).^2./(2*sig.^2))./sqrt(2*pi*sig.^2) + q/2*exp(-(d18O(n)-mu-d*sig).^2./(2*sig.^2))./sqrt(2*pi*sig.^2) + q/2*exp(-(d18O(n)-mu+d*sig).^2./(2*sig.^2))./sqrt(2*pi*sig.^2);
                % MargLik = MargLik + log(AA);
                MargLik = MargLik - (a_d18O+0.5)*log(1+(d18O(n)-mu).^2./(2*b_d18O*sig.^2)) - log(sig) + log_tail;
            end
        end
        if strcmp(mode,'d18O') == 0 && ~isempty(C14_Table)
            MargLik = MargLik + EP_Y(SAM_A(3,:),C14_Table,cal_curve,a,b);
        end
        
        if ~isnan(Age_Info(1))
            if Age_Info(3) == 0
                MargLik = MargLik - (SAM_A(3,:)-Age_Info(1)).^2./(2*Age_Info(2).^2);
            elseif Age_Info(3) == 1
                index = (SAM_A(3,:)>=Age_Info(1)-Age_Info(2))&(SAM_A(3,:)<=Age_Info(1)+Age_Info(2));
                MargLik(~index) = -inf;
            end
        end
        
        MargLik(SAM_A(3,:)<age_st) = -inf;
        MargLik(isnan(MargLik)) = -inf;
        
        if ~isnan(data.max)
            MargLik(SAM_A(3,:)>data.max) = -inf;
        end
        if ~isnan(data.min)
            MargLik(SAM_A(3,:)<data.min) = -inf;
        end
        
        WW(3,:) = MargLik;
    else
        index = (~isinf(WWB(1,:)));
        AB_TABLE_C = AB(1,index);
        WB_TABLE_C = WWB(1,index);
        
        index = (~isinf(WWB(2,:)));
        AB_TABLE_M = AB(2,index);
        WB_TABLE_M = WWB(2,index);
        
        index = (~isinf(WWB(3,:)));
        AB_TABLE_E = AB(3,index);
        WB_TABLE_E = WWB(3,index);
        
        % AB = AB(:,1:min(S,200));
        % WWB = WWB(:,1:min(S,200));
        %{
        AB_TABLE_C = AB(1,:);
        AB_TABLE_M = AB(2,:);
        AB_TABLE_E = AB(3,:);
        
        WB_TABLE_C = WWB(1,:);
        WB_TABLE_M = WWB(2,:);
        WB_TABLE_E = WWB(3,:);
        %}
        
        if ~isempty(data.initial_age)
            ID = (abs(data.initial_age(:,1)-data.depth(nn))<1e-6);
            if sum(ID) > 0
                ZZ = data.initial_age(ID,2);
                if ZZ(1) > age_st && ZZ(1) < age_ed
                    SAM_A = ZZ(1)*ones(3,S);
                else
                    SAM_A = (age_ed-age_st).*rand(3,S) + age_st;
                end
            else
                % Sample ages:
                SAM_A = (age_ed-age_st).*rand(3,S) + age_st;
            end
        else
            % Sample ages:
            SAM_A = (age_ed-age_st).*rand(3,S) + age_st;
        end
        
        % Compute weights:
        WW = zeros(3,S);
        
        RR_C = interp1(R(:,1),R(:,2),AB_TABLE_C');
        RR_C = repmat(RR_C,[1,S]);
        
        RR_M = interp1(R(:,1),R(:,2),AB_TABLE_M');
        RR_M = repmat(RR_M,[1,S]);
        
        RR_E = interp1(R(:,1),R(:,2),AB_TABLE_E');
        RR_E = repmat(RR_E,[1,S]);
        
        % Z == 1:
        MargLik_C = PHI(1,1) + WB_TABLE_C';
        MargLik_M = PHI(2,1) + WB_TABLE_M';
        MargLik_E = PHI(3,1) + WB_TABLE_E';
        
        VV_C = (AB_TABLE_C'-SAM_A(1,:))./(RR_C.*depth_diff);
        VV_M = (AB_TABLE_M'-SAM_A(1,:))./(RR_M.*depth_diff);
        VV_E = (AB_TABLE_E'-SAM_A(1,:))./(RR_E.*depth_diff);
        
        MargLik_C = MargLik_C + interp1(data.ACC_MODEL(:,1),data.ACC_MODEL(:,2),VV_C,'linear',-56) - data.ACC_CONTRACTION;
        index = (VV_C<=0)|(VV_C>1./1.0850);
        MargLik_C(index) = -inf;
        
        MargLik_M = MargLik_M + interp1(data.ACC_MODEL(:,1),data.ACC_MODEL(:,2),VV_M,'linear',-56) - data.ACC_CONTRACTION;
        index = (VV_M<=0)|(VV_M>1./1.0850);
        MargLik_M(index) = -inf;
        
        MargLik_E = MargLik_E + interp1(data.ACC_MODEL(:,1),data.ACC_MODEL(:,2),VV_E,'linear',-56) - data.ACC_CONTRACTION;
        index = (VV_E<=0)|(VV_E>1./1.0850);
        MargLik_E(index) = -inf;
        
        index = (VV_C>1./data.lower_sedrate+0.001)|(VV_C<1./data.upper_sedrate-0.001);
        MargLik_C(index) = -inf;
        
        index = (VV_M>1./data.lower_sedrate+0.001)|(VV_M<1./data.upper_sedrate-0.001);
        MargLik_M(index) = -inf;
        
        index = (VV_E>1./data.lower_sedrate+0.001)|(VV_E<1./data.upper_sedrate-0.001);
        MargLik_E(index) = -inf;
        
        MargLik = [MargLik_C;MargLik_M;MargLik_E];
        
        AMAX = max(MargLik);
        MargLik = AMAX + log(sum(exp(MargLik-AMAX)));
        
        if strcmp(mode,'C14') == 0 && ~isnan(d18O(1))
            % MargLik = MargLik + log((1-q)*exp(EP_V(SAM_A(1,:),zeros(1,S),d18O,stack))+q*exp(EP_V(SAM_A(1,:),ones(1,S),d18O,stack)));
            mu = interp1(stack(:,1),stack(:,2),SAM_A(1,:));
            sig = interp1(stack(:,1),stack(:,3),SAM_A(1,:));
            N = sum(~isnan(d18O));
            for n = 1:N
                % AA = (1-q)*exp(-(d18O(n)-mu).^2./(2*sig.^2))./sqrt(2*pi*sig.^2) + q/2*exp(-(d18O(n)-mu-d*sig).^2./(2*sig.^2))./sqrt(2*pi*sig.^2) + q/2*exp(-(d18O(n)-mu+d*sig).^2./(2*sig.^2))./sqrt(2*pi*sig.^2);
                % MargLik = MargLik + log(AA);
                MargLik = MargLik - (a_d18O+0.5)*log(1+(d18O(n)-mu).^2./(2*b_d18O*sig.^2)) - log(sig) + log_tail;
            end
        end
        if strcmp(mode,'d18O') == 0 && ~isempty(C14_Table)
            MargLik = MargLik + EP_Y(SAM_A(1,:),C14_Table,cal_curve,a,b);
        end
        
        if ~isnan(Age_Info(1))
            if Age_Info(3) == 0
                MargLik = MargLik - (SAM_A(1,:)-Age_Info(1)).^2./(2*Age_Info(2).^2);
            elseif Age_Info(3) == 1
                index = (SAM_A(1,:)>=Age_Info(1)-Age_Info(2))&(SAM_A(1,:)<=Age_Info(1)+Age_Info(2));
                MargLik(~index) = -inf;
            end
        end
        
        MargLik(SAM_A(1,:)<age_st) = -inf;
        MargLik(isnan(MargLik)) = -inf;
        
        if ~isnan(data.max)
            MargLik(SAM_A(1,:)>data.max) = -inf;
        end
        if ~isnan(data.min)
            MargLik(SAM_A(1,:)<data.min) = -inf;
        end
        
        WW(1,:) = MargLik;
        
        % Z == 2:
        MargLik_C = PHI(1,2) + WB_TABLE_C';
        MargLik_M = PHI(2,2) + WB_TABLE_M';
        MargLik_E = PHI(3,2) + WB_TABLE_E';
        
        VV_C = (AB_TABLE_C'-SAM_A(2,:))./(RR_C.*depth_diff);
        VV_M = (AB_TABLE_M'-SAM_A(2,:))./(RR_M.*depth_diff);
        VV_E = (AB_TABLE_E'-SAM_A(2,:))./(RR_E.*depth_diff);
        
        MargLik_C = MargLik_C + interp1(data.ACC_MODEL(:,1),data.ACC_MODEL(:,2),VV_C,'linear',-56) - data.ACC_STEADY;
        index = (VV_C<=1./1.0850)|(VV_C>1./0.9220);
        MargLik_C(index) = -inf;
        
        MargLik_M = MargLik_M + interp1(data.ACC_MODEL(:,1),data.ACC_MODEL(:,2),VV_M,'linear',-56) - data.ACC_STEADY;
        index = (VV_M<=1./1.0850)|(VV_M>1./0.9220);
        MargLik_M(index) = -inf;
        
        MargLik_E = MargLik_E + interp1(data.ACC_MODEL(:,1),data.ACC_MODEL(:,2),VV_E,'linear',-56) - data.ACC_STEADY;
        index = (VV_E<=1./1.0850)|(VV_E>1./0.9220);
        MargLik_E(index) = -inf;
        
        index = (VV_C>1./data.lower_sedrate+0.001)|(VV_C<1./data.upper_sedrate-0.001);
        MargLik_C(index) = -inf;
        
        index = (VV_M>1./data.lower_sedrate+0.001)|(VV_M<1./data.upper_sedrate-0.001);
        MargLik_M(index) = -inf;
        
        index = (VV_E>1./data.lower_sedrate+0.001)|(VV_E<1./data.upper_sedrate-0.001);
        MargLik_E(index) = -inf;
        
        MargLik = [MargLik_C;MargLik_M;MargLik_E];
        
        AMAX = max(MargLik);
        MargLik = AMAX + log(sum(exp(MargLik-AMAX)));
        
        if strcmp(mode,'C14') == 0 && ~isnan(d18O(1))
            % MargLik = MargLik + log((1-q)*exp(EP_V(SAM_A(2,:),zeros(1,S),d18O,stack))+q*exp(EP_V(SAM_A(2,:),ones(1,S),d18O,stack)));
            mu = interp1(stack(:,1),stack(:,2),SAM_A(2,:));
            sig = interp1(stack(:,1),stack(:,3),SAM_A(2,:));
            N = sum(~isnan(d18O));
            for n = 1:N
                % AA = (1-q)*exp(-(d18O(n)-mu).^2./(2*sig.^2))./sqrt(2*pi*sig.^2) + q/2*exp(-(d18O(n)-mu-d*sig).^2./(2*sig.^2))./sqrt(2*pi*sig.^2) + q/2*exp(-(d18O(n)-mu+d*sig).^2./(2*sig.^2))./sqrt(2*pi*sig.^2);
                % MargLik = MargLik + log(AA);
                MargLik = MargLik - (a_d18O+0.5)*log(1+(d18O(n)-mu).^2./(2*b_d18O*sig.^2)) - log(sig) + log_tail;
            end
        end
        if strcmp(mode,'d18O') == 0 && ~isempty(C14_Table)
            MargLik = MargLik + EP_Y(SAM_A(2,:),C14_Table,cal_curve,a,b);
        end
        
        if ~isnan(Age_Info(1))
            if Age_Info(3) == 0
                MargLik = MargLik - (SAM_A(2,:)-Age_Info(1)).^2./(2*Age_Info(2).^2);
            elseif Age_Info(3) == 1
                index = (SAM_A(2,:)>=Age_Info(1)-Age_Info(2))&(SAM_A(2,:)<=Age_Info(1)+Age_Info(2));
                MargLik(~index) = -inf;
            end
        end
        
        MargLik(SAM_A(2,:)<age_st) = -inf;
        MargLik(isnan(MargLik)) = -inf;
        
        if ~isnan(data.max)
            MargLik(SAM_A(2,:)>data.max) = -inf;
        end
        if ~isnan(data.min)
            MargLik(SAM_A(2,:)<data.min) = -inf;
        end
        
        WW(2,:) = MargLik;
        
        % Z == 3:
        MargLik_C = PHI(1,3) + WB_TABLE_C';
        MargLik_M = PHI(2,3) + WB_TABLE_M';
        MargLik_E = PHI(3,3) + WB_TABLE_E';
        
        VV_C = (AB_TABLE_C'-SAM_A(3,:))./(RR_C.*depth_diff);
        VV_M = (AB_TABLE_M'-SAM_A(3,:))./(RR_M.*depth_diff);
        VV_E = (AB_TABLE_E'-SAM_A(3,:))./(RR_E.*depth_diff);
        
        MargLik_C = MargLik_C + interp1(data.ACC_MODEL(:,1),data.ACC_MODEL(:,2),VV_C,'linear',-56) - data.ACC_EXPANSION;
        index = (VV_C<=1./0.9220);
        MargLik_C(index) = -inf;
        
        MargLik_M = MargLik_M + interp1(data.ACC_MODEL(:,1),data.ACC_MODEL(:,2),VV_M,'linear',-56) - data.ACC_EXPANSION;
        index = (VV_M<1./0.9220);
        MargLik_M(index) = -inf;
        
        MargLik_E = MargLik_E + interp1(data.ACC_MODEL(:,1),data.ACC_MODEL(:,2),VV_E,'linear',-56) - data.ACC_EXPANSION;
        index = (VV_E<1./0.9220);
        MargLik_E(index) = -inf;
        
        index = (VV_C>1./data.lower_sedrate+0.001)|(VV_C<1./data.upper_sedrate-0.001);
        MargLik_C(index) = -inf;
        
        index = (VV_M>1./data.lower_sedrate+0.001)|(VV_M<1./data.upper_sedrate-0.001);
        MargLik_M(index) = -inf;
        
        index = (VV_E>1./data.lower_sedrate+0.001)|(VV_E<1./data.upper_sedrate-0.001);
        MargLik_E(index) = -inf;
        
        MargLik = [MargLik_C;MargLik_M;MargLik_E];
        
        AMAX = max(MargLik);
        MargLik = AMAX + log(sum(exp(MargLik-AMAX)));
        
        if strcmp(mode,'C14') == 0 && ~isnan(d18O(1))
            % MargLik = MargLik + log((1-q)*exp(EP_V(SAM_A(3,:),zeros(1,S),d18O,stack))+q*exp(EP_V(SAM_A(3,:),ones(1,S),d18O,stack)));
            mu = interp1(stack(:,1),stack(:,2),SAM_A(3,:));
            sig = interp1(stack(:,1),stack(:,3),SAM_A(3,:));
            N = sum(~isnan(d18O));
            for n = 1:N
                % AA = (1-q)*exp(-(d18O(n)-mu).^2./(2*sig.^2))./sqrt(2*pi*sig.^2) + q/2*exp(-(d18O(n)-mu-d*sig).^2./(2*sig.^2))./sqrt(2*pi*sig.^2) + q/2*exp(-(d18O(n)-mu+d*sig).^2./(2*sig.^2))./sqrt(2*pi*sig.^2);
                % MargLik = MargLik + log(AA);
                MargLik = MargLik - (a_d18O+0.5)*log(1+(d18O(n)-mu).^2./(2*b_d18O*sig.^2)) - log(sig) + log_tail;
            end
        end
        if strcmp(mode,'d18O') == 0 && ~isempty(C14_Table)
            MargLik = MargLik + EP_Y(SAM_A(3,:),C14_Table,cal_curve,a,b);
        end
        
        if ~isnan(Age_Info(1))
            if Age_Info(3) == 0
                MargLik = MargLik - (SAM_A(3,:)-Age_Info(1)).^2./(2*Age_Info(2).^2);
            elseif Age_Info(3) == 1
                index = (SAM_A(3,:)>=Age_Info(1)-Age_Info(2))&(SAM_A(3,:)<=Age_Info(1)+Age_Info(2));
                MargLik(~index) = -inf;
            end
        end
        
        MargLik(SAM_A(3,:)<age_st) = -inf;
        MargLik(isnan(MargLik)) = -inf;
        
        if ~isnan(data.max)
            MargLik(SAM_A(3,:)>data.max) = -inf;
        end
        if ~isnan(data.min)
            MargLik(SAM_A(3,:)<data.min) = -inf;
        end
        
        WW(3,:) = MargLik;
    end
end


end