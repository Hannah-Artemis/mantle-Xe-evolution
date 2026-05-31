function[data_5,data_25,data_50,data_75,data_95] = ...
calculate_percentile_fun(data,nt,t,size)
% this function returns the middle 50%, 90%, and the median
% values of the given dataset
%
% modified on Meng Guo, Jun Korenaga 
% create empty matrixes to store the results
data_5 = zeros(1,nt);
data_25 = zeros(1,nt);
data_50 = zeros(1,nt);
data_75 = zeros(1,nt);
data_95 = zeros(1,nt);
% sort the given dataset
for i= 1:nt
    data_sort(i,:) = sort(data(i,:));
end
% get the middle 50%, 90%, and median values of the given dataset
for i = 1:nt
    data_25(i) = data_sort(i,max(1,round(size*0.25)));
    data_75(i) = data_sort(i,max(1,round(size*0.75)));
    data_5(i) = data_sort(i,max(1,round(size*0.05)));
    data_95(i) = data_sort(i,max(1,round(size*0.95)));
    data_50(i) = data_sort(i,max(1,round(size*0.50)));
end