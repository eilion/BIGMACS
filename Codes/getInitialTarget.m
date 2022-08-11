function [target] = getInitialTarget(inputFile,setting,MODE)

target = struct('stack',cell(1,1),'cal_curve',cell(1,1));


% load stack:
stack_path = ['Inputs/',inputFile,'/stack.txt'];
if exist(stack_path,'file') ~= 2
    stack_path = 'Defaults/stack.txt';
end

fileID = fopen(stack_path);
INFO = textscan(fileID,'%s %s %s');
fclose(fileID);

if isempty(INFO{2}) == 1
    path = ['Defaults/Regional_Stacks/',INFO{1}{1},'.txt'];
    fileID = fopen(path);
    STACK = textscan(fileID,'%s %s %s');
    fclose(fileID);
    
    target.stack = zeros(length(STACK{1})-1,3);
    
    AGE = STACK{1}(2:end);
    MEAN = STACK{2}(2:end);
    STDV = STACK{3}(2:end);
    
    target.stack(:,1) = str2double(AGE);
    target.stack(:,2) = str2double(MEAN);
    target.stack(:,3) = str2double(STDV);
else
    target.stack = zeros(length(INFO{1})-1,3);
    
    AGE = INFO{1}(2:end);
    MEAN = INFO{2}(2:end);
    STDV = INFO{3}(2:end);
    
    target.stack(:,1) = str2double(AGE);
    target.stack(:,2) = str2double(MEAN);
    target.stack(:,3) = str2double(STDV);
end

INDEX = (target.stack(:,1)>=setting.stack_min)&(target.stack(:,1)<=setting.stack_max);
target.stack = target.stack(INDEX,:);

if strcmp(MODE,'stacking')
    st = max([setting.st,target.stack(1,1)]);
    ed = min([setting.ed,target.stack(end,1)]);
    age = (st:setting.interval:ed)';
    
    index = (age>=target.stack(1,1))&(age<=target.stack(end,1));
    age = age(index);
    
    N = length(age);
    new_stack = zeros(N,3);
    new_stack(:,1) = age;
    new_stack(:,2) = interp1(target.stack(:,1),target.stack(:,2),age);
    new_stack(:,3) = interp1(target.stack(:,1),target.stack(:,3),age);
    
    target.stack = new_stack;
end

target.init_stack = target.stack;


% load calibration curve:
cal_curve = cell(3,1);
cal_curve{1} = load('Defaults/Calibration_Curves/IntCal20.txt');
cal_curve{2} = load('Defaults/Calibration_Curves/Marine20.txt');
cal_curve{3} = load('Defaults/Calibration_Curves/SHCal20.txt');

cal_curve_ms = cell(3,1);
age = (0:0.05:target.stack(end,1))';

if age(end) < 60
    age = [age;60];
end

for m = 1:3
    
    Cal_age = cal_curve{m}(:,1)/1000;
    C14_age = cal_curve{m}(:,2)/1000;
    C14_error = cal_curve{m}(:,5)/1000;
    
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

target.cal_curve = cal_curve_ms;


end