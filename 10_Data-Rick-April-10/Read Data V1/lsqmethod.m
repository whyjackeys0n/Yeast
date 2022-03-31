clear
clc

load CENPK21C_0410.mat

CENPK21C_diff_sucrose_pb = CENPK21C_0410.CENPK21C_diff_sucrose_pb;
T = CENPK21C_diff_sucrose_pb(1,4:end);


ytrue = CENPK21C_diff_sucrose_pb(11,4:end);

f = @(c,x)c(1)./(1+exp(-c(2)*(T-2*c(3))));

c = lsqcurvefit(f,[1e5 1 5],T,ytrue);

plot(T,ytrue,'o-',T,f(c,T),'r.:')