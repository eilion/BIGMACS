function [stack] = getStack(stack_path)

fileID = fopen(stack_path);
INFO = textscan(fileID,'%s %s %s');
fclose(fileID);

if isempty(INFO{2}) == 1
    path = ['Defaults/Regional_Stacks/',INFO{1}{1},'.txt'];
    fileID = fopen(path);
    STACK = textscan(fileID,'%s %s %s');
    fclose(fileID);
    
    stack = zeros(length(STACK{1})-1,3);
    
    AGE = STACK{1}(2:end);
    MEAN = STACK{2}(2:end);
    STDV = STACK{3}(2:end);
    
    stack(:,1) = str2double(AGE);
    stack(:,2) = str2double(MEAN);
    stack(:,3) = str2double(STDV);
else
    stack = zeros(length(INFO{1})-1,3);
    
    AGE = INFO{1}(2:end);
    MEAN = INFO{2}(2:end);
    STDV = INFO{3}(2:end);
    
    stack(:,1) = str2double(AGE);
    stack(:,2) = str2double(MEAN);
    stack(:,3) = str2double(STDV);
end


end

