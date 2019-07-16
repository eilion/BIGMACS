addpath('Codes/');

% compare two age inferences:
titles = cell(2,1);
results_path = cell(2,1);
titles{1} = 'Dual Proxy Ages';
titles{2} = '14C Ages';
results_path{1} = 'Outputs/Stack_Construction/DNEA_stack_both_OU';
results_path{2} = 'Outputs/Core_Alignments/DNEA_stack_C14/';

AgeVsAge(titles,results_path);