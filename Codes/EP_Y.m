function [log_prob] = EP_Y(A,C14_Table,cal_curve,a,b)

log_tail = log(gamma(a+0.5)) - log(gamma(a)) - 0.5*log(2*pi*b);

T = size(A,2);
L = size(C14_Table,1);

log_prob = zeros(1,T);

if isempty(C14_Table) == 0
    for ll = 1:L
        c14_mu = interp1(cal_curve{C14_Table(ll,5)}(:,1),cal_curve{C14_Table(ll,5)}(:,2),A);
        c14_stdv = interp1(cal_curve{C14_Table(ll,5)}(:,1),cal_curve{C14_Table(ll,5)}(:,3),A);
        
        if C14_Table(ll,end) == 1
            log_prob = log_prob - 0.5*(c14_mu+C14_Table(ll,3)-C14_Table(ll,1)).^2./(c14_stdv.^2+C14_Table(ll,2)^2+C14_Table(ll,4)^2) - 0.5*log(c14_stdv.^2+C14_Table(ll,2)^2+C14_Table(ll,4)^2) - 0.5*log(2*pi);
        else
            log_prob = log_prob - (a+0.5)*log(1+(c14_mu+C14_Table(ll,3)-C14_Table(ll,1)).^2./(2*b*(c14_stdv.^2+C14_Table(ll,2)^2+C14_Table(ll,4)^2))) - 0.5*log(c14_stdv.^2+C14_Table(ll,2)^2+C14_Table(ll,4)^2) + log_tail;
        end
    end
end

log_prob(isnan(log_prob)) = 0;


end