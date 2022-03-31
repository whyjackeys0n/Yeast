clear
clc

load CENPK21C_0410.mat

CENPK21C_diff_sucrose_pb = CENPK21C_0410.CENPK21C_diff_sucrose_pb;
T = CENPK21C_diff_sucrose_pb(1,4:end);
r = [];
K = [];
L = [];

for i = 1:3
    ytrue = CENPK21C_diff_sucrose_pb(i+1,4:end);
    f = @(c,x)c(1)./(1+exp(-c(2)*(T-c(3))));
    c = lsqcurvefit(f,[1e5 1 10],T,ytrue);
    K = [K,c(1)];
    r = [r,c(2)];
    L = [L,c(3)];
    plot(T,ytrue);
    f1 = plot(T,f(c,T),'r.');
    hold on
end

for i = 4:6
    ytrue = CENPK21C_diff_sucrose_pb(i+1,4:end);
    f = @(c,x)c(1)./(1+exp(-c(2)*(T-c(3))));
    c = lsqcurvefit(f,[1e5 1 10],T,ytrue);
    K = [K,c(1)];
    r = [r,c(2)];
    L = [L,c(3)];
    plot(T,ytrue);
    f2 = plot(T,f(c,T),'b.');
    hold on
end

for i = 7:9
    ytrue = CENPK21C_diff_sucrose_pb(i+1,4:end);
    f = @(c,x)c(1)./(1+exp(-c(2)*(T-c(3))));
    c = lsqcurvefit(f,[1e5 1 10],T,ytrue);
    K = [K,c(1)];
    r = [r,c(2)];
    L = [L,c(3)];
    plot(T,ytrue);
    f3 = plot(T,f(c,T),'c.');
    hold on
end

for i = 10:12
    ytrue = CENPK21C_diff_sucrose_pb(i+1,4:end);
    f = @(c,x)c(1)./(1+exp(-c(2)*(T-c(3))));
    c = lsqcurvefit(f,[1e5 1 5],T,ytrue);
    K = [K,c(1)];
    r = [r,c(2)];
    L = [L,c(3)];
    plot(T,ytrue);
    f4 = plot(T,f(c,T),'g.');
    hold on
end

legend([f1,f2,f3,f4],'Sucrose Concentration: 116.9 mM','Sucrose Concentration: 29.2 mM','Sucrose Concentration: 7.3 mM','Sucrose Concentration: 1.83 mM','FontSize',18,'Location','SouthEast');

set(gcf,'position',[300,300,1500,1000]);

