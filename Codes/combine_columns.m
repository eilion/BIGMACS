function column = combine_columns(depth, d18O_data)

% Combines d18O data to one column
% Input 1: d18O depth
% Input 2: d18O data matrix
% Output: [depth, d18O]
% 8/08/2019

% # of rows
sz1 = length(d18O_data(:,1));
% # of columns
sz2 = length(d18O_data(1,:));

% loops through # of columns to combine to one
for i = 1:sz2
    column((sz1*(i-1))+1:sz1*i,1) = depth;
    column((sz1*(i-1))+1:sz1*i,2) = d18O_data(:,i);
end


% deletes NaN rows
column(isnan(column(:,2))==1,:)=[];
% sort by depth
column = sortrows(column,1);
    