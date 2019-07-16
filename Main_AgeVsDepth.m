addpath('Codes/');

% plot single alignments:
results_path = cell(1,1);
results_path{1} = 'Outputs/Stack_Construction/EXAMPLE_INPUT_both_OU';

%{
% compare multiple alignments:
results_path = cell(2,1);
results_path{1} = 'Outputs/Core_Alignments/Two_Cores_DNEA_both';
results_path{2} = 'Outputs/Core_Alignments/Two_Cores_C14';
%}

AgeVsDepth(results_path);