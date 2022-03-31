clear;
clc;

load CENPK21C_0410.mat
load S288C_0410.mat

S288C_diff_sucrose_pbprch = S288C_0410.S288C_diff_sucrose_pbprch;

time = S288C_diff_sucrose_pbprch(1,124:end)'-28;
c = S288C_diff_sucrose_pbprch(2,124:end)';




x0 = [20000, 1, 2];

f = fittype('K/(1+exp(-r*(t-L)))','independent','t','coefficients',{'K','r','L'}); 
cfun=fit(time,c,f);
xi=0:0.1:100;

plot(time,c,'r*');
hold on
[X,resnorm]=lsqcurvefit(@GrowthRate,x0,time,c);
y2i = X(1)./(1+exp(-X(2).*(xi-X(3))));
plot(xi,y2i)