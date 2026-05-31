function[data_low,data_up] = ...
calculate_margin_fun(data,nt,t,size)
% this function returns the middle 50%, 90%, and the median
% values of the given dataset
%
% modified on (Guo & Korenaga, 2020) 
% create empty matrixes to store the results
data_low = zeros(1,nt);
data_up = zeros(1,nt);

% sort the given dataset
for i= 1:nt
    data_sort(i,:) = sort(data(i,:));
end
% get the middle 50%, 90%, and median values of the given dataset
for i = 1:nt
    data_low(i) = data_sort(i,max(1,round(size*0)));
    data_up(i) = data_sort(i,max(1,round(size*1)));
end