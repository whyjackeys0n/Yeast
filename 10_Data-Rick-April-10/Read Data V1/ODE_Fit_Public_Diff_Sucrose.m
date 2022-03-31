clear
clc

load CENPK21C_0410.mat

CENPK21C_diff_sucrose_pb = CENPK21C_0410.CENPK21C_diff_sucrose_pb;
T = CENPK21C_diff_sucrose_pb(1,4:end)';
r = [];
K = [];
for i = 1:12
    ytrue = CENPK21C_diff_sucrose_pb(i+1,4:end)';

    x0 = [15 1 1e5];   
    ub = [1000 10 1e7];        
    lb = [0 0 0];         
    F = @(x) COST(x,T,ytrue);
    xout = fmincon(F,x0,[],[],[],[],lb,ub); 
    r = [r,xout(2)];
    K = [K,xout(3)];
end

function COST = COST(x,T,ytrue)

y0 = x(1);
r = x(2);
K = x(3);
% The cost function calls the ODE solver.
[tout,yout] = ode45(@dydt,T,y0,[],r,K);
COST = sum((yout - ytrue).^2);

end

function yprime = dydt(t,y,r,K)
yprime = r*y*(1-y/K);
end