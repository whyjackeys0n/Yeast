clear
clc

load CENPK21C_0410.mat
load S288C_0410.mat

S288C_diff_sucrose_pbprch = S288C_0410.S288C_diff_sucrose_pbprch;
T = S288C_diff_sucrose_pbprch(1,2:end)';
ytrue = S288C_diff_sucrose_pbprch(5,2:end)';

plot(T,ytrue);
hold on;

h = plot(T,nan*ytrue,'r');
set(h,'tag','solution');
x0 = [15 1 1e5];   
ub = [100 10 1e6];        
lb = [0 0 0];         
F = @(x) COST(x,T,ytrue);
xout = fmincon(F,x0,[],[],[],[],lb,ub); 
legend({'Experimental Data','Fitted Data'});
xout(2)


function COST = COST(x,T,ytrue)


y0 = x(1);
r = x(2);
K = x(3);
% The cost function calls the ODE solver.
[tout,yout] = ode45(@dydt,T,y0,[],r,K);
COST = sum((yout - ytrue).^2);
h = findobj('tag','solution');
set(h,'ydata',yout);
title(['y0 = ' num2str(y0) '   ' 'r = ' num2str(r) '    K = ' num2str(K)]);

drawnow;


end

function yprime = dydt(t,y,r,K)
yprime = r*y*(1-y/K);
end
