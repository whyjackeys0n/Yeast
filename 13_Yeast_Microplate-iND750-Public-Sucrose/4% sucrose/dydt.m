function yprime = dydt(t,y,r,K)
yprime = r*y*(1-y/K);
end