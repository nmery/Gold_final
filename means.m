function y = means(x)
%
% MEANS  Average or mean value. For column vectors MEANS(x) returns the mean
%        value. For matrices or row vector, MEANS(x) is a row vector containing
%        the mean value of each column. The only difference with MATLAB
%        function MEAN is for a row vector, where MEANS returns the row vector
%        instead of the mean value of the elements of the row.
%
[m,n] = size(x);
if m > 1,
   y = sum(x) / m;
else
   y = x;
end
