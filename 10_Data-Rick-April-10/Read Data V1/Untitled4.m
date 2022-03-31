function xout = do_optimization
% The true parameters
y0_true = 1;
A_true = 2;
B_true = 3;
T = (0:.01:1)';
ytrue =3./(1 + 2*exp(-6*T)) + 0.1*randn(size(T));
hold on;
plot(T,ytrue);
h = plot(T,nan*ytrue,'r');
set(h,'tag','solution');
x0 = [0.4 3.9 1.2];   % Just some Initial Condition
ub = [5 5 5];         % Upper bounds
lb = [0 0 0];         % Lower bounds
F = @(x) COST(x,T,ytrue);
xout = fmincon(F,x0,[],[],[],[],lb,ub); %<-- FMINCON is the optimizer
legend({'Experimental Data','Fitted Data'});

function COST = COST(x,T,ytrue)
y0 = x(1);
A = x(2);
B = x(3);
% The cost function calls the ODE solver.
[tout,yout] = ode45(@dydt,T,y0,[],A,B);
COST = sum((yout - ytrue).^2);
h = findobj('tag','solution');
set(h,'ydata',yout);
title(['y0 = ' num2str(y0) '   ' 'A = ' num2str(A) '    B = ' num2str(B)]);
drawnow;

function yprime = dydt(t,y,A,B)
yprime = A*y*(B-y);