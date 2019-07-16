function [SAM_A,WW] = Proposal_5(WWB,AB,depth_diff,d18O,C14_Table,Age_Info,core_param,global_param,stack,cal_curve,mode,QQ)

alpha = global_param.alpha;
beta = global_param.beta;

S = global_param.nParticles;
v = global_param.v;

R = core_param.R;

phi_I = core_param.phi_I;
phi_C = core_param.phi_C;
phi_M = core_param.phi_M;
phi_E = core_param.phi_E;
PHI = [phi_C;phi_M;phi_E;phi_I];
PHI = log(PHI);

q = global_param.q;
d = global_param.d;

if isempty(WWB) == 1
    % Sample ages:
    if size(QQ,2) == S
        index = ceil(S*rand(1,5*S));
        QQ = QQ(index);
        UU = rand(1,5*S)*(4*v) + (QQ-2*v);
        UU = UU((UU>=core_param.min)&(UU<=core_param.max));
        SAM_A = UU(1:S);
    else
        UU = rand(1,5*S)*(4*QQ(2)) + (QQ(1)-2*QQ(2));
        UU = UU((UU>=core_param.min)&(UU<=core_param.max));
        SAM_A = UU(1:S);
    end
    
    % Compute weights:
    MargLik = zeros(1,S);
    if strcmp(mode,'C14') == 0 && isnan(d18O(1)) == 0
        mu = interp1(stack(:,1),stack(:,2),SAM_A);
        sig = interp1(stack(:,1),stack(:,3),SAM_A);
        N = sum(isnan(d18O)==0);
        for n = 1:N
            AA = (1-q)*exp(-(d18O(n)-mu).^2./(2*sig.^2))./sqrt(2*pi*sig.^2) + q/2*exp(-(d18O(n)-mu-d*sig).^2./(2*sig.^2))./sqrt(2*pi*sig.^2) + q/2*exp(-(d18O(n)-mu+d*sig).^2./(2*sig.^2))./sqrt(2*pi*sig.^2);
            MargLik = MargLik + log(AA);
        end
    end
    if strcmp(mode,'d18O') == 0
        MargLik = MargLik + EP_Y(SAM_A,C14_Table,cal_curve);
    end
    
    if isnan(Age_Info(1)) == 0
        if Age_Info(3) == 0
            MargLik = MargLik - (SAM_A-Age_Info(1)).^2./(2*Age_Info(2).^2);
        elseif Age_Info(3) == 1
            index = (SAM_A>Age_Info(1)-Age_Info(2))&(SAM_A<Age_Info(1)+Age_Info(2));
            MargLik(~index) = -inf;
        end
    end
    
    MargLik(SAM_A<core_param.min|SAM_A>core_param.max) = -inf;
    MargLik(isnan(MargLik)==1) = -inf;
    
    if isnan(core_param.max) == 0
        MargLik(SAM_A>core_param.max) = -inf;
    end
    if isnan(core_param.min) == 0
        MargLik(SAM_A<core_param.min) = -inf;
    end
    
    WW = MargLik;
else
    if size(WWB,1) == 1
        index = (isinf(WWB)==0);
        AB_TABLE = AB(index);
        WB_TABLE = WWB(index);
        
        % AB_TABLE = AB(1:min(asamples,S));
        % WB_TABLE = WWB(1:min(asamples,S));
        
        % Sample ages:
        if size(QQ,2) == S
            index = ceil(S*rand(1,15*S));
            QQ = QQ(index);
            UU = rand(1,15*S)*(4*v) + (QQ-2*v);
            UU = UU((UU>=core_param.min)&(UU<=core_param.max));
            SAM_A = reshape(UU(1:3*S),[3,S]);
        else
            UU = rand(1,15*S)*(4*QQ(2)) + (QQ(1)-2*QQ(2));
            UU = UU((UU>=core_param.min)&(UU<=core_param.max));
            SAM_A = reshape(UU(1:3*S),[3,S]);
        end
        
        % Compute weights:
        WW = zeros(3,S);
        
        RR = interp1(R(:,1),R(:,2),AB_TABLE');
        RR = repmat(RR,[1,S]);
        
        % Z == 1:
        MargLik = PHI(4,1) + WB_TABLE';
        
        VV = (AB_TABLE'-SAM_A(1,:))./(RR.*depth_diff);
        
        MargLik = MargLik + (alpha-1)*log(VV) - beta*VV - log(RR) - log(gamcdf(0.9220,alpha,1/beta));
        index = (VV<=0)|(VV>=0.9220);
        MargLik(index) = -inf;
        
        index = (VV<core_param.lower_sedrate)|(VV>core_param.upper_sedrate);
        MargLik(index) = -inf;
        
        AMAX = max(MargLik);
        MargLik = AMAX + log(sum(exp(MargLik-AMAX)));
        
        if strcmp(mode,'C14') == 0 && isnan(d18O(1)) == 0
            % MargLik = MargLik + log((1-q)*exp(EP_V(SAM_A(1,:),zeros(1,S),d18O,stack))+q*exp(EP_V(SAM_A(1,:),ones(1,S),d18O,stack)));
            mu = interp1(stack(:,1),stack(:,2),SAM_A(1,:));
            sig = interp1(stack(:,1),stack(:,3),SAM_A(1,:));
            N = sum(isnan(d18O)==0);
            for n = 1:N
                AA = (1-q)*exp(-(d18O(n)-mu).^2./(2*sig.^2))./sqrt(2*pi*sig.^2) + q/2*exp(-(d18O(n)-mu-d*sig).^2./(2*sig.^2))./sqrt(2*pi*sig.^2) + q/2*exp(-(d18O(n)-mu+d*sig).^2./(2*sig.^2))./sqrt(2*pi*sig.^2);
                MargLik = MargLik + log(AA);
            end
        end
        if strcmp(mode,'d18O') == 0
            MargLik = MargLik + EP_Y(SAM_A(1,:),C14_Table,cal_curve);
        end
        
        if isnan(Age_Info(1)) == 0
            if Age_Info(3) == 0
                MargLik = MargLik - (SAM_A(1,:)-Age_Info(1)).^2./(2*Age_Info(2).^2);
            elseif Age_Info(3) == 1
                index = (SAM_A(1,:)>Age_Info(1)-Age_Info(2))&(SAM_A(1,:)<Age_Info(1)+Age_Info(2));
                MargLik(~index) = -inf;
            end
        end
        %{
        if size(QQ,2) == S
            MargLik = MargLik + (a+0.5)*log(1+(SAM_A(1,:)-QQ(1,:)).^2./(2*b*v^2));
        else
            MargLik = MargLik + (a+0.5)*log(1+(SAM_A(1,:)-QQ(1)).^2./(2*b*QQ(2)));
        end
        %}
        MargLik(SAM_A(1,:)<core_param.min) = -inf;
        MargLik(isnan(MargLik)==1) = -inf;
        
        if isnan(core_param.max) == 0
            MargLik(SAM_A(1,:)>core_param.max) = -inf;
        end
        if isnan(core_param.min) == 0
            MargLik(SAM_A(1,:)<core_param.min) = -inf;
        end
        
        WW(1,:) = MargLik;
        
        % Z == 2:
        MargLik = PHI(4,2) + WB_TABLE';
        
        VV = (AB_TABLE'-SAM_A(2,:))./(RR.*depth_diff);
        
        MargLik = MargLik + (alpha-1)*log(VV) - beta*VV - log(RR) - log(gamcdf(1.0850,alpha,1/beta)-gamcdf(0.9220,alpha,1/beta));
        index = (VV<0.9220)|(VV>=1.0850);
        MargLik(index) = -inf;
        
        index = (VV<core_param.lower_sedrate)|(VV>core_param.upper_sedrate);
        MargLik(index) = -inf;
        
        AMAX = max(MargLik);
        MargLik = AMAX + log(sum(exp(MargLik-AMAX)));
        
        if strcmp(mode,'C14') == 0 && isnan(d18O(1)) == 0
            % MargLik = MargLik + log((1-q)*exp(EP_V(SAM_A(2,:),zeros(1,S),d18O,stack))+q*exp(EP_V(SAM_A(2,:),ones(1,S),d18O,stack)));
            mu = interp1(stack(:,1),stack(:,2),SAM_A(2,:));
            sig = interp1(stack(:,1),stack(:,3),SAM_A(2,:));
            N = sum(isnan(d18O)==0);
            for n = 1:N
                AA = (1-q)*exp(-(d18O(n)-mu).^2./(2*sig.^2))./sqrt(2*pi*sig.^2) + q/2*exp(-(d18O(n)-mu-d*sig).^2./(2*sig.^2))./sqrt(2*pi*sig.^2) + q/2*exp(-(d18O(n)-mu+d*sig).^2./(2*sig.^2))./sqrt(2*pi*sig.^2);
                MargLik = MargLik + log(AA);
            end
        end
        if strcmp(mode,'d18O') == 0
            MargLik = MargLik + EP_Y(SAM_A(2,:),C14_Table,cal_curve);
        end
        
        if isnan(Age_Info(1)) == 0
            if Age_Info(3) == 0
                MargLik = MargLik - (SAM_A(2,:)-Age_Info(1)).^2./(2*Age_Info(2).^2);
            elseif Age_Info(3) == 1
                index = (SAM_A(2,:)>Age_Info(1)-Age_Info(2))&(SAM_A(2,:)<Age_Info(1)+Age_Info(2));
                MargLik(~index) = -inf;
            end
        end
        %{
        if size(QQ,2) == S
            MargLik = MargLik + (a+0.5)*log(1+(SAM_A(2,:)-QQ(2,:)).^2./(2*b*v^2));
        else
            MargLik = MargLik + (a+0.5)*log(1+(SAM_A(2,:)-QQ(1)).^2./(2*b*QQ(2)));
        end
        %}
        MargLik(SAM_A(2,:)<core_param.min) = -inf;
        MargLik(isnan(MargLik)==1) = -inf;
        
        if isnan(core_param.max) == 0
            MargLik(SAM_A(2,:)>core_param.max) = -inf;
        end
        if isnan(core_param.min) == 0
            MargLik(SAM_A(2,:)<core_param.min) = -inf;
        end
        
        WW(2,:) = MargLik;
        
        % Z == 3:
        MargLik = PHI(4,3) + WB_TABLE';
        
        VV = (AB_TABLE'-SAM_A(3,:))./(RR.*depth_diff);
        
        MargLik = MargLik + (alpha-1)*log(VV) - beta*VV - log(RR) - log(1-gamcdf(1.0850,alpha,1/beta));
        index = (VV<1.0850);
        MargLik(index) = -inf;
        
        index = (VV<core_param.lower_sedrate)|(VV>core_param.upper_sedrate);
        MargLik(index) = -inf;
        
        AMAX = max(MargLik);
        MargLik = AMAX + log(sum(exp(MargLik-AMAX)));
        
        if strcmp(mode,'C14') == 0 && isnan(d18O(1)) == 0
            % MargLik = MargLik + log((1-q)*exp(EP_V(SAM_A(3,:),zeros(1,S),d18O,stack))+q*exp(EP_V(SAM_A(3,:),ones(1,S),d18O,stack)));
            mu = interp1(stack(:,1),stack(:,2),SAM_A(3,:));
            sig = interp1(stack(:,1),stack(:,3),SAM_A(3,:));
            N = sum(isnan(d18O)==0);
            for n = 1:N
                AA = (1-q)*exp(-(d18O(n)-mu).^2./(2*sig.^2))./sqrt(2*pi*sig.^2) + q/2*exp(-(d18O(n)-mu-d*sig).^2./(2*sig.^2))./sqrt(2*pi*sig.^2) + q/2*exp(-(d18O(n)-mu+d*sig).^2./(2*sig.^2))./sqrt(2*pi*sig.^2);
                MargLik = MargLik + log(AA);
            end
        end
        if strcmp(mode,'d18O') == 0
            MargLik = MargLik + EP_Y(SAM_A(3,:),C14_Table,cal_curve);
        end
        
        if isnan(Age_Info(1)) == 0
            if Age_Info(3) == 0
                MargLik = MargLik - (SAM_A(3,:)-Age_Info(1)).^2./(2*Age_Info(2).^2);
            elseif Age_Info(3) == 1
                index = (SAM_A(3,:)>Age_Info(1)-Age_Info(2))&(SAM_A(3,:)<Age_Info(1)+Age_Info(2));
                MargLik(~index) = -inf;
            end
        end
        %{
        if size(QQ,2) == S
            MargLik = MargLik + (a+0.5)*log(1+(SAM_A(3,:)-QQ(1,:)).^2./(2*b*v^2));
        else
            MargLik = MargLik + (a+0.5)*log(1+(SAM_A(3,:)-QQ(1)).^2./(2*b*QQ(2)));
        end
        %}
        MargLik(SAM_A(3,:)<core_param.min) = -inf;
        MargLik(isnan(MargLik)==1) = -inf;
        
        if isnan(core_param.max) == 0
            MargLik(SAM_A(3,:)>core_param.max) = -inf;
        end
        if isnan(core_param.min) == 0
            MargLik(SAM_A(3,:)<core_param.min) = -inf;
        end
        
        WW(3,:) = MargLik;
    else
        index = (isinf(WWB(1,:))==0);
        AB_TABLE_C = AB(1,index);
        WB_TABLE_C = WWB(1,index);
        
        index = (isinf(WWB(2,:))==0);
        AB_TABLE_M = AB(2,index);
        WB_TABLE_M = WWB(2,index);
        
        index = (isinf(WWB(3,:))==0);
        AB_TABLE_E = AB(3,index);
        WB_TABLE_E = WWB(3,index);
        
        % Sample ages:
        if size(QQ,2) == S
            index = ceil(S*rand(1,15*S));
            QQ = QQ(index);
            UU = rand(1,15*S)*(4*v) + (QQ-2*v);
            UU = UU((UU>=core_param.min)&(UU<=core_param.max));
            SAM_A = reshape(UU(1:3*S),[3,S]);
        else
            UU = rand(1,15*S)*(4*QQ(2)) + (QQ(1)-2*QQ(2));
            UU = UU((UU>=core_param.min)&(UU<=core_param.max));
            SAM_A = reshape(UU(1:3*S),[3,S]);
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
        
        MargLik_C = MargLik_C + (alpha-1)*log(VV_C) - beta*VV_C - log(RR_C) - log(gamcdf(0.9220,alpha,1/beta));
        index = (VV_C<=0)|(VV_C>=0.9220);
        MargLik_C(index) = -inf;
        
        MargLik_M = MargLik_M + (alpha-1)*log(VV_M) - beta*VV_M - log(RR_M) - log(gamcdf(0.9220,alpha,1/beta));
        index = (VV_M<=0)|(VV_M>=0.9220);
        MargLik_M(index) = -inf;
        
        MargLik_E = MargLik_E + (alpha-1)*log(VV_E) - beta*VV_E - log(RR_E) - log(gamcdf(0.9220,alpha,1/beta));
        index = (VV_E<=0)|(VV_E>=0.9220);
        MargLik_E(index) = -inf;
        
        index = (VV_C<core_param.lower_sedrate)|(VV_C>core_param.upper_sedrate);
        MargLik_C(index) = -inf;
        
        index = (VV_M<core_param.lower_sedrate)|(VV_M>core_param.upper_sedrate);
        MargLik_M(index) = -inf;
        
        index = (VV_E<core_param.lower_sedrate)|(VV_E>core_param.upper_sedrate);
        MargLik_E(index) = -inf;
        
        MargLik = [MargLik_C;MargLik_M;MargLik_E];
        
        AMAX = max(MargLik);
        MargLik = AMAX + log(sum(exp(MargLik-AMAX)));
        
        if strcmp(mode,'C14') == 0 && isnan(d18O(1)) == 0
            % MargLik = MargLik + log((1-q)*exp(EP_V(SAM_A(1,:),zeros(1,S),d18O,stack))+q*exp(EP_V(SAM_A(1,:),ones(1,S),d18O,stack)));
            mu = interp1(stack(:,1),stack(:,2),SAM_A(1,:));
            sig = interp1(stack(:,1),stack(:,3),SAM_A(1,:));
            N = sum(isnan(d18O)==0);
            for n = 1:N
                AA = (1-q)*exp(-(d18O(n)-mu).^2./(2*sig.^2))./sqrt(2*pi*sig.^2) + q/2*exp(-(d18O(n)-mu-d*sig).^2./(2*sig.^2))./sqrt(2*pi*sig.^2) + q/2*exp(-(d18O(n)-mu+d*sig).^2./(2*sig.^2))./sqrt(2*pi*sig.^2);
                MargLik = MargLik + log(AA);
            end
        end
        if strcmp(mode,'d18O') == 0
            MargLik = MargLik + EP_Y(SAM_A(1,:),C14_Table,cal_curve);
        end
        
        if isnan(Age_Info(1)) == 0
            if Age_Info(3) == 0
                MargLik = MargLik - (SAM_A(1,:)-Age_Info(1)).^2./(2*Age_Info(2).^2);
            elseif Age_Info(3) == 1
                index = (SAM_A(1,:)>Age_Info(1)-Age_Info(2))&(SAM_A(1,:)<Age_Info(1)+Age_Info(2));
                MargLik(~index) = -inf;
            end
        end
        %{
        if size(QQ,2) == S
            MargLik = MargLik + (a+0.5)*log(1+(SAM_A(1,:)-QQ(1,:)).^2./(2*b*v^2));
        else
            MargLik = MargLik + (a+0.5)*log(1+(SAM_A(1,:)-QQ(1)).^2./(2*b*QQ(2)));
        end
        %}
        MargLik(SAM_A(1,:)<core_param.min) = -inf;
        MargLik(isnan(MargLik)==1) = -inf;
        
        if isnan(core_param.max) == 0
            MargLik(SAM_A(1,:)>core_param.max) = -inf;
        end
        if isnan(core_param.min) == 0
            MargLik(SAM_A(1,:)<core_param.min) = -inf;
        end
        
        WW(1,:) = MargLik;
        
        % Z == 2:
        MargLik_C = PHI(1,2) + WB_TABLE_C';
        MargLik_M = PHI(2,2) + WB_TABLE_M';
        MargLik_E = PHI(3,2) + WB_TABLE_E';
        
        VV_C = (AB_TABLE_C'-SAM_A(2,:))./(RR_C.*depth_diff);
        VV_M = (AB_TABLE_M'-SAM_A(2,:))./(RR_M.*depth_diff);
        VV_E = (AB_TABLE_E'-SAM_A(2,:))./(RR_E.*depth_diff);
        
        MargLik_C = MargLik_C + (alpha-1)*log(VV_C) - beta*VV_C - log(RR_C) - log(gamcdf(1.0850,alpha,1/beta)-gamcdf(0.9220,alpha,1/beta));
        index = (VV_C<0.9220)|(VV_C>=1.0850);
        MargLik_C(index) = -inf;
        
        MargLik_M = MargLik_M + (alpha-1)*log(VV_M) - beta*VV_M - log(RR_M) - log(gamcdf(1.0850,alpha,1/beta)-gamcdf(0.9220,alpha,1/beta));
        index = (VV_M<0.9220)|(VV_M>=1.0850);
        MargLik_M(index) = -inf;
        
        MargLik_E = MargLik_E + (alpha-1)*log(VV_E) - beta*VV_E - log(RR_E) - log(gamcdf(1.0850,alpha,1/beta)-gamcdf(0.9220,alpha,1/beta));
        index = (VV_E<0.9220)|(VV_E>=1.0850);
        MargLik_E(index) = -inf;
        
        index = (VV_C<core_param.lower_sedrate)|(VV_C>core_param.upper_sedrate);
        MargLik_C(index) = -inf;
        
        index = (VV_M<core_param.lower_sedrate)|(VV_M>core_param.upper_sedrate);
        MargLik_M(index) = -inf;
        
        index = (VV_E<core_param.lower_sedrate)|(VV_E>core_param.upper_sedrate);
        MargLik_E(index) = -inf;
        
        MargLik = [MargLik_C;MargLik_M;MargLik_E];
        
        AMAX = max(MargLik);
        MargLik = AMAX + log(sum(exp(MargLik-AMAX)));
        
        if strcmp(mode,'C14') == 0 && isnan(d18O(1)) == 0
            % MargLik = MargLik + log((1-q)*exp(EP_V(SAM_A(2,:),zeros(1,S),d18O,stack))+q*exp(EP_V(SAM_A(2,:),ones(1,S),d18O,stack)));
            mu = interp1(stack(:,1),stack(:,2),SAM_A(2,:));
            sig = interp1(stack(:,1),stack(:,3),SAM_A(2,:));
            N = sum(isnan(d18O)==0);
            for n = 1:N
                AA = (1-q)*exp(-(d18O(n)-mu).^2./(2*sig.^2))./sqrt(2*pi*sig.^2) + q/2*exp(-(d18O(n)-mu-d*sig).^2./(2*sig.^2))./sqrt(2*pi*sig.^2) + q/2*exp(-(d18O(n)-mu+d*sig).^2./(2*sig.^2))./sqrt(2*pi*sig.^2);
                MargLik = MargLik + log(AA);
            end
        end
        if strcmp(mode,'d18O') == 0
            MargLik = MargLik + EP_Y(SAM_A(2,:),C14_Table,cal_curve);
        end
        
        if isnan(Age_Info(1)) == 0
            if Age_Info(3) == 0
                MargLik = MargLik - (SAM_A(2,:)-Age_Info(1)).^2./(2*Age_Info(2).^2);
            elseif Age_Info(3) == 1
                index = (SAM_A(2,:)>Age_Info(1)-Age_Info(2))&(SAM_A(2,:)<Age_Info(1)+Age_Info(2));
                MargLik(~index) = -inf;
            end
        end
        %{
        if size(QQ,2) == S
            MargLik = MargLik + (a+0.5)*log(1+(SAM_A(2,:)-QQ(2,:)).^2./(2*b*v^2));
        else
            MargLik = MargLik + (a+0.5)*log(1+(SAM_A(2,:)-QQ(1)).^2./(2*b*QQ(2)));
        end
        %}
        MargLik(SAM_A(2,:)<core_param.min) = -inf;
        MargLik(isnan(MargLik)==1) = -inf;
        
        if isnan(core_param.max) == 0
            MargLik(SAM_A(2,:)>core_param.max) = -inf;
        end
        if isnan(core_param.min) == 0
            MargLik(SAM_A(2,:)<core_param.min) = -inf;
        end
        
        WW(2,:) = MargLik;
        
        % Z == 3:
        MargLik_C = PHI(1,3) + WB_TABLE_C';
        MargLik_M = PHI(2,3) + WB_TABLE_M';
        MargLik_E = PHI(3,3) + WB_TABLE_E';
        
        VV_C = (AB_TABLE_C'-SAM_A(3,:))./(RR_C.*depth_diff);
        VV_M = (AB_TABLE_M'-SAM_A(3,:))./(RR_M.*depth_diff);
        VV_E = (AB_TABLE_E'-SAM_A(3,:))./(RR_E.*depth_diff);
        
        MargLik_C = MargLik_C + (alpha-1)*log(VV_C) - beta*VV_C - log(RR_C) - log(1-gamcdf(1.0850,alpha,1/beta));
        index = (VV_C<1.0850);
        MargLik_C(index) = -inf;
        
        MargLik_M = MargLik_M + (alpha-1)*log(VV_M) - beta*VV_M - log(RR_M) - log(1-gamcdf(1.0850,alpha,1/beta));
        index = (VV_M<1.0850);
        MargLik_M(index) = -inf;
        
        MargLik_E = MargLik_E + (alpha-1)*log(VV_E) - beta*VV_E - log(RR_E) - log(1-gamcdf(1.0850,alpha,1/beta));
        index = (VV_E<1.0850);
        MargLik_E(index) = -inf;
        
        index = (VV_C<core_param.lower_sedrate)|(VV_C>core_param.upper_sedrate);
        MargLik_C(index) = -inf;
        
        index = (VV_M<core_param.lower_sedrate)|(VV_M>core_param.upper_sedrate);
        MargLik_M(index) = -inf;
        
        index = (VV_E<core_param.lower_sedrate)|(VV_E>core_param.upper_sedrate);
        MargLik_E(index) = -inf;
        
        MargLik = [MargLik_C;MargLik_M;MargLik_E];
        
        AMAX = max(MargLik);
        MargLik = AMAX + log(sum(exp(MargLik-AMAX)));
        
        if strcmp(mode,'C14') == 0 && isnan(d18O(1)) == 0
            % MargLik = MargLik + log((1-q)*exp(EP_V(SAM_A(3,:),zeros(1,S),d18O,stack))+q*exp(EP_V(SAM_A(3,:),ones(1,S),d18O,stack)));
            mu = interp1(stack(:,1),stack(:,2),SAM_A(3,:));
            sig = interp1(stack(:,1),stack(:,3),SAM_A(3,:));
            N = sum(isnan(d18O)==0);
            for n = 1:N
                AA = (1-q)*exp(-(d18O(n)-mu).^2./(2*sig.^2))./sqrt(2*pi*sig.^2) + q/2*exp(-(d18O(n)-mu-d*sig).^2./(2*sig.^2))./sqrt(2*pi*sig.^2) + q/2*exp(-(d18O(n)-mu+d*sig).^2./(2*sig.^2))./sqrt(2*pi*sig.^2);
                MargLik = MargLik + log(AA);
            end
        end
        if strcmp(mode,'d18O') == 0
            MargLik = MargLik + EP_Y(SAM_A(3,:),C14_Table,cal_curve);
        end
        
        if isnan(Age_Info(1)) == 0
            if Age_Info(3) == 0
                MargLik = MargLik - (SAM_A(3,:)-Age_Info(1)).^2./(2*Age_Info(2).^2);
            elseif Age_Info(3) == 1
                index = (SAM_A(3,:)>Age_Info(1)-Age_Info(2))&(SAM_A(3,:)<Age_Info(1)+Age_Info(2));
                MargLik(~index) = -inf;
            end
        end
        %{
        if size(QQ,2) == S
            MargLik = MargLik + (a+0.5)*log(1+(SAM_A(3,:)-QQ(3,:)).^2./(2*b*v^2));
        else
            MargLik = MargLik + (a+0.5)*log(1+(SAM_A(3,:)-QQ(1)).^2./(2*b*QQ(2)));
        end
        %}
        MargLik(SAM_A(3,:)<core_param.min) = -inf;
        MargLik(isnan(MargLik)==1) = -inf;
        
        if isnan(core_param.max) == 0
            MargLik(SAM_A(3,:)>core_param.max) = -inf;
        end
        if isnan(core_param.min) == 0
            MargLik(SAM_A(3,:)<core_param.min) = -inf;
        end
        
        WW(3,:) = MargLik;
    end
end


end