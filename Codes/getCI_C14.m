function [CI] = getCI_C14(data,target)

a = 3;
b = 4;

cal_curve = target.cal_curve;

L = length(data);
CI = struct('name',cell(L,1),'depth',cell(L,1),'lower',cell(L,1),'median',cell(L,1),'upper',cell(L,1),'samples',cell(L,1));

parfor ll = 1:L
    CI(ll).name = data(ll).name;
    M = 0;
    for n = 1:length(data(ll).depth)
        if isempty(data(ll).radiocarbon{n}) == 0
            M = M + size(data(ll).radiocarbon{n},1);
        end
    end
    
    depth = zeros(M,1);
    lower = zeros(M,2);
    MEDIAN = zeros(M,2);
    upper = zeros(M,2);
    SAMPLES = zeros(M,1000);
    
    k = 0;
    for n = 1:length(data(ll).depth)
        if isempty(data(ll).radiocarbon{n}) == 0
            Table = data(ll).radiocarbon{n};
            for m = 1:size(Table,1)
                k = k + 1;
                depth(k) = data(ll).depth(n);
                
                rand_seed = log(rand(100,10000));
                
                c14_age = min(max(Table(m,1)-Table(m,3),cal_curve{Table(m,5)}(1,2)),cal_curve{Table(m,5)}(end,2));
                x = interp1q(cal_curve{Table(m,5)}(:,2),cal_curve{Table(m,5)}(:,1),c14_age)*ones(1,10000);
                
                for iters = 1:100
                    c14_x = min(max(x,cal_curve{Table(m,5)}(1,1)),cal_curve{Table(m,5)}(end,1));
                    det_x = interp1q(cal_curve{Table(m,5)}(:,1),cal_curve{Table(m,5)}(:,2),c14_x')';
                    err_x = interp1q(cal_curve{Table(m,5)}(:,1),cal_curve{Table(m,5)}(:,3),c14_x')';
                    RR_x = - (a+0.5)*log(2*b+(det_x+Table(m,3)-Table(m,1)).^2./(err_x.^2+Table(m,2)^2+Table(m,4)^2)) - 0.5*log(err_x.^2+Table(m,2)^2+Table(m,4)^2);
                    
                    z = normrnd(x,0.3);
                    c14_z = min(max(z,cal_curve{Table(m,5)}(1,1)),cal_curve{Table(m,5)}(end,1));                    
                    det_z = interp1q(cal_curve{Table(m,5)}(:,1),cal_curve{Table(m,5)}(:,2),c14_z')';
                    err_z = interp1q(cal_curve{Table(m,5)}(:,1),cal_curve{Table(m,5)}(:,3),c14_z')';
                    RR_z = - (a+0.5)*log(2*b+(det_z+Table(m,3)-Table(m,1)).^2./(err_z.^2+Table(m,2)^2+Table(m,4)^2)) - 0.5*log(err_z.^2+Table(m,2)^2+Table(m,4)^2);
                    
                    RR_z(z<0|z>55) = -inf;
                    
                    index = (RR_z-RR_x>rand_seed(iters,:));
                    x(index) = z(index);
                end
                
                x_samples = x;
                SAMPLES(k,:) = x_samples(1:1000);
                
                lower(k,1) = quantile(x_samples,0.025);
                lower(k,2) = quantile(x_samples,0.16);
                MEDIAN(k,1) = quantile(x_samples,0.5);
                MEDIAN(k,2) = mean(x_samples)
                upper(k,1) = quantile(x_samples,0.84);
                upper(k,2) = quantile(x_samples,0.975);
            end
        end
    end
    CI(ll).depth = depth;
    CI(ll).lower = lower;
    CI(ll).median = MEDIAN;
    CI(ll).upper = upper;
    CI(ll).samples = SAMPLES;
end


end