function [param,target] = standardizeStack(data,param,target)

L = length(data);

DET = 0;
for ll = 1:L
    if strcmp(data(ll).islearn_shift,'no') && strcmp(data(ll).islearn_scale,'no')
        DET = 1;
    end
end

if DET == 0
    for ll = 1:L
        if strcmp(data(ll).islearn_shift,'no')
            DET = 2;
        end
    end
end

if DET == 0
    for ll = 1:L
        if strcmp(data(ll).islearn_scale,'no')
            DET = 3;
        end
    end
end

%{
if DET == 1
    stack(:,2) = stack(:,2)*global_param.scale + global_param.shift;
    stack(:,3) = stack(:,3)*global_param.scale;
    
    STACK_SAMPLE(:,:,1) = STACK_SAMPLE(:,:,1)*global_param.scale + global_param.shift;
    STACK_SAMPLE(:,:,2) = STACK_SAMPLE(:,:,2)*global_param.scale + global_param.shift;
    
    for ll = 1:L
        rcd_stack(ll).mean = stack(:,2);
        rcd_stack(ll).stdv = stack(:,3);
    end
elseif DET == 2
    stack(:,2) = stack(:,2) + global_param.shift/global_param.scale;
    
    STACK_SAMPLE(:,:,1) = STACK_SAMPLE(:,:,1) + global_param.shift/global_param.scale;
    STACK_SAMPLE(:,:,2) = STACK_SAMPLE(:,:,2) + global_param.shift/global_param.scale;
    
    for ll = 1:L
        rcd_stack(ll).mean = stack(:,2);
    end
    
    for ll = 1:L
        core_param(ll).scale = core_param(ll).scale*global_param.scale;
    end
elseif DET == 3
    stack(:,2) = stack(:,2)*global_param.scale;
    stack(:,3) = stack(:,3)*global_param.scale;
    
    STACK_SAMPLE(:,:,1) = STACK_SAMPLE(:,:,1)*global_param.scale;
    STACK_SAMPLE(:,:,2) = STACK_SAMPLE(:,:,2)*global_param.scale;
    
    for ll = 1:L
        rcd_stack(ll).mean = stack(:,2);
        rcd_stack(ll).stdv = stack(:,3);
    end
    
    for ll = 1:L
        core_param(ll).shift = core_param(ll).scale*global_param.shift + core_param(ll).shift;
    end
else
    for ll = 1:L
        core_param(ll).shift = core_param(ll).scale*global_param.shift + core_param(ll).shift;
        core_param(ll).scale = core_param(ll).scale*global_param.scale;
    end
end
%}

if DET == 1 || DET == 2 || DET == 3
    target.stack(:,2) = target.stack(:,2)*param.scale + param.shift;
    target.stack(:,3) = target.stack(:,3)*param.scale;
    
    target.stack_sample(:,:,1) = target.stack_sample(:,:,1)*param.scale + param.shift;
    target.stack_sample(:,:,2) = target.stack_sample(:,:,2)*param.scale + param.shift;
else
    AA = 0;
    BB = 0;
    for ll = 1:L
        AA = AA + data(ll).scale;
        BB = BB + data(ll).shift;
    end
    AA = AA/L;
    BB = BB/L;
    param.scale = param.scale*AA;
    param.shift = param.shift*AA + BB;
    for ll = 1:L
        data(ll).shift = data(ll).shift - data(ll).scale*BB/AA;
        data(ll).scale = data(ll).scale/AA;
    end
    
    target.stack(:,2) = target.stack(:,2)*param.scale + param.shift;
    target.stack(:,3) = target.stack(:,3)*param.scale;
    
    target.stack_sample(:,:,1) = target.stack_sample(:,:,1)*param.scale + param.shift;
    target.stack_sample(:,:,2) = target.stack_sample(:,:,2)*param.scale + param.shift;
end

param.scale = 1;
param.shift = 0;


end