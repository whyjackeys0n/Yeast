%% main program
clear
clc

load CENPK21C_0410.mat

CENPK21C_diff_sucrose_pb = CENPK21C_0410.CENPK21C_diff_sucrose_pb;
T = CENPK21C_diff_sucrose_pb(1,4:end)';
ytrue = CENPK21C_diff_sucrose_pb(2,4:end)';
figure(1)
set(gcf,'position',[300,300,1500,1000]);
plot(T,ytrue);
hold on;

h = plot(T,nan*ytrue,'r');
set(h,'tag','solution');
x0 = [15 1 1e5];   
ub = [1000 10 1e7];        
lb = [0 0 0];         
F = @(x) COST(x,T,ytrue);
xout = fmincon(F,x0,[],[],[],[],lb,ub); 
legend({'Experimental Data','Fitted Data'},'FontSize',20,'Location','SouthEast');
xout(2)

%% cost function
function COST = COST(x,T,ytrue)
y0 = x(1);
r = x(2);
K = x(3);
% The cost function calls the ODE solver.
[tout,yout] = ode45(@dydt,T,y0,[],r,K);
COST = sum((yout - ytrue).^2);
h = findobj('tag','solution');
set(h,'ydata',yout);
title(['y0 = ' num2str(y0) '   ' 'r = ' num2str(r) '    K = ' num2str(K)],'FontSize',20);
set(gcf,'position',[300,300,1500,1000]);
drawnow;
end

%% ordinary differential equation
function yprime = dydt(t,y,r,K)
yprime = r*y*(1-y/K);
end
