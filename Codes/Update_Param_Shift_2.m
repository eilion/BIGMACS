function [core_param] = Update_Param_Shift_2(data,Samples,core_param,stack,mode)

beta1 = 0.9;
beta2 = 0.999;
epsilon = 1e-8;
gamma = 1e-4;

if strcmp(mode,'d18O') == 1 || strcmp(mode,'both') == 1
    L = length(data);
    
    AGE = stack(:,1);
    MU = stack(:,2);
    STDV = stack(:,3);
    
    parfor ll = 1:L
        
        [~,~,~,islearn_shift,islearn_scale,~,~,~] = getSetting_Core(data(ll).name);
        
        if strcmp(islearn_shift,'yes') == 1 || strcmp(islearn_scale,'yes') == 1
            
            M = size(Samples(ll).ages,2);
            V = cell(M,1);
            A = cell(M,1);
            
            for m = 1:M
                index = Samples(ll).isoutlier{m};
                V{m} = data(ll).d18O(index==0);
                A{m} = repmat(Samples(ll).ages(:,m),[1,size(data(ll).d18O,2)]);
                A{m} = A{m}(index==0);
            end
            
            V = cat(1,V{:});
            A = cat(1,A{:});
            
            mu = interp1(AGE,MU,A);
            sig = interp1(AGE,STDV,A);
            
            h = core_param(ll).shift;
            c = core_param(ll).scale;
            Mw = zeros(2,1);
            Vw = zeros(2,1);
            PDEV = zeros(2,1);
            for r = 1:10000
                
                if strcmp(islearn_shift,'yes') == 1
                    PDEV(1) = sum((V-c*mu-h).*(c*sig).^(-2));
                end
                if strcmp(islearn_scale,'yes') == 1
                    PDEV(2) = sum((V-c*mu-h).*(c*sig).^(-3).*(V-h).*sig-c.^(-1));
                end
                
                Mw = beta1*Mw - (1-beta1).*PDEV;
                Vw = beta2*Vw + (1-beta2).*PDEV.*PDEV;
                
                h = h - (gamma*sqrt(1-beta2^r)/(1-beta1^r)).*Mw(1)./(sqrt(Vw(1))+epsilon);
                c = c - (gamma*sqrt(1-beta2^r)/(1-beta1^r)).*Mw(2)./(sqrt(Vw(2))+epsilon);
            end
            
            core_param(ll).shift = h;
            core_param(ll).scale = c;
        end
    end
    
end


end

