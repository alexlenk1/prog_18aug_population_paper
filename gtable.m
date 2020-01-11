
function y = table(x,ndec)

% x is a matrix containing the data to be displayed
% ndec is the number of decimal places to display (with rounding)
% returns a matrix of characters formatted for a LaTeX table
%
% Keisuke Hirano
% khirano@fas.harvard.edu
% 1 May 1998

[n,k] = size(x);
entrystring = strcat(' %1.',num2str(round(ndec)));
entrystring = strcat(entrystring,'f');
formatstring = ' & ';
formatstring = strcat(formatstring,entrystring);
for j=2:k,
    formatstring = strcat(formatstring,' & ');
    formatstring = strcat(formatstring,entrystring);
end;
formatstring = strcat(formatstring,' \\\\');

for i=1:n,
    % using two steps seems easier since different lines can have 
    % different lengths
    tempstring = sprintf(formatstring,x(i,:));
    y(i,1:length(tempstring)) = tempstring;
end;

% end table.m

