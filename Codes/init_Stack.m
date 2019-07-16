function [new_stack] = init_Stack(stack,age)

index = (age>=stack(1,1))&(age<=stack(end,1));
age = age(index);

N = length(age);
new_stack = zeros(N,3);
new_stack(:,1) = age;
new_stack(:,2) = interp1(stack(:,1),stack(:,2),age);
new_stack(:,3) = interp1(stack(:,1),stack(:,3),age);


end

