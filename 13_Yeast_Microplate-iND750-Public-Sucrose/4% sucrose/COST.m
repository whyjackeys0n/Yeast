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

