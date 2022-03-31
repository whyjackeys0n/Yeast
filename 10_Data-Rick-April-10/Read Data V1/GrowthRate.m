function dN = GrowthRate(param,t)
dN = param(1)./(1+exp(-param(2).*(t-param(3))));