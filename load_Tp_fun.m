function[Tp_anchorHerz1,Tp_anchorHerz2,Tp_anchorHerz3,Tp_anchorHerz4,...
t_anchorHerz1,t_anchorHerz2,t_anchorHerz3,t_anchorHerz4,...
t_Herz,Tp_Herz] = load_Tp_fun(data_Tp,t_pd)

% this function load in the anchor points for Herzburg's data
%
% modified on (Guo & Korenaga, 2020) 
% column 1 is the age backward in time, column 2 is the mantle ...
% potential temperature,column 3 is the age forward in time with original
% tmax
t_Herz = t_pd-data_Tp(:,1);% age forward in time 
Tp_Herz = data_Tp(:,2);
% chose four anchor points of Herzburg's data to make a trend
Tp_anchorHerz1 = 1470; t_Herz1back=0.8;
Tp_anchorHerz2 = 1554; t_Herz2back=1.87;
Tp_anchorHerz3 = 1578; t_Herz3back=2.75;
Tp_anchorHerz4 = 1567; t_Herz4back=3.39;
t_anchorHerz1 = t_pd-t_Herz1back;
t_anchorHerz2 = t_pd-t_Herz2back;
t_anchorHerz3 = t_pd-t_Herz3back;
t_anchorHerz4 = t_pd-t_Herz4back;