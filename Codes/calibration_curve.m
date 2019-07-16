function [cal_curve_ms] = calibration_curve(age)

cal_curve = cell(3,1);
cal_curve{1} = load('Calibration_Curves/IntCal13.txt');
cal_curve{2} = load('Calibration_Curves/Marine13.txt');
cal_curve{3} = load('Calibration_Curves/SHCal13.txt');

cal_curve_ms = cell(3,1);

if age(end) < 60
    age = [age;60];
end

for m = 1:3
    
    Cal_age = cal_curve{m}(:,1)/1000;
    C14_age = cal_curve{m}(:,2)/1000;
    C14_error = cal_curve{m}(:,5)/1000;
    
    % Cal_age_add = [Cal_age;Cal_age(end)+10];
    % C14_age_add = [C14_age;10*(C14_age(end)-C14_age(1))/(Cal_age(end)-Cal_age(1))+C14_age(end)];
    % C14_error_add = [C14_error;10*(C14_error(end)-C14_error(1))/(Cal_age(end)-Cal_age(1))+C14_error(end)];
    Cal_age_add = [Cal_age;age(end)];
    C14_age_add = [C14_age;(age(end)-Cal_age(end))*(C14_age(end)-C14_age(1))/(Cal_age(end)-Cal_age(1))+C14_age(end)];
    C14_error_add = [C14_error;(age(end)-Cal_age(end))*(C14_error(end)-C14_error(1))/(Cal_age(end)-Cal_age(1))+C14_error(end)];
    
    cal_curve_ms{m} = zeros(length(age),3);
    cal_curve_ms{m}(:,1) = age;
    
    index = (cal_curve_ms{m}(:,1) >= Cal_age_add(1) & cal_curve_ms{m}(:,1) <= Cal_age_add(end));
    cal_curve_ms{m}(index,2) = interp1(Cal_age_add,C14_age_add,cal_curve_ms{m}(index,1));
    cal_curve_ms{m}(index,3) = interp1(Cal_age_add,C14_error_add,cal_curve_ms{m}(index,1));
    
    cal_curve_ms{m}(~index,2) = NaN;
    cal_curve_ms{m}(~index,3) = NaN;
end


end

