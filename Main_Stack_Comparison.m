addpath('Codes/');

% compare multiple stacks:
results_path = cell(4,1);
results_path{1} = 'Outputs/Stack_Construction/DNEA_stack_both_OU';
results_path{2} = 'Outputs/Stack_Construction/DNEA_stack_both_M15';
results_path{3} = 'Outputs/Stack_Construction/DNEA_stack_both_M25';
results_path{4} = 'Outputs/Stack_Construction/DNEA_stack_both_SE';

titles = cell(4,1);
titles{1} = 'DNEA_Stack (OU)';
titles{2} = 'DNEA_Stack (M15)';
titles{3} = 'DNEA_Stack (M25)';
titles{4} = 'DNEA_Stack (SE)';


Stack_Comparison(titles,results_path);