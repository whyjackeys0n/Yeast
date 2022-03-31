function [lb,ub] = RHS(t,y,INFO)

% This subroutine updates the upper and lower bounds for the fluxes in the
% exID arrays in main. The output should be two matrices, lb and ub. The lb matrix
% contains the lower bounds for exID{i} in the ith row in the same order as
% exID. The same is true for the upper bounds in the ub matrix.
% Infinity can be used for unconstrained variables, however, it should be 
% fixed for all time. 

%% Assign value
param = INFO.param;

VmaxG   = param(1);
KmG     = param(2);
VmaxF   = param(3);
KmF     = param(4);
VmaxS   = param(5);
KmS     = param(6);
VmaxO   = param(7);
KmO     = param(8);
VmaxPho = param(9);       
KmPho   = param(10);            
VmaxSul = param(11);          
KmSul   = param(12);

%% Set extracellular metabolite concentrations
glc = y(4);
fru = y(5);
scr = y(6);
eth = y(7);
oxy = y(9);
pho = y(10);
sul = y(11);

% Public - can uptake glucose and fructose but not sucrose
lb(1,1) = 0;
ub(1,1) = Inf;

lb(1,2) = -VmaxG * max([glc, 0]) / (KmG + glc);
ub(1,2) = Inf;

lb(1,3) = -VmaxF * max([fru, 0]) / (KmF + fru);
ub(1,3) = Inf;

lb(1,4) = 0;
ub(1,4) = Inf;

lb(1,5) = 0;
ub(1,5) = Inf;

lb(1,6) = -VmaxO * max([oxy, 0]) / (KmO + oxy);
ub(1,6) = Inf;

lb(1,7) = 0;
ub(1,7) = Inf;

lb(1,8) = 0;
ub(1,8) = Inf;

% Private – can uptake glucose, fructose and sucrose
lb(2,1) = 0;
ub(2,1) = Inf;

lb(2,2) = -VmaxG * max([glc, 0]) / (KmG + glc);
ub(2,2) = Inf;

lb(2,3) = -VmaxF * max([fru, 0]) / (KmF+ fru);
ub(2,3) = Inf;

%lb(2,4) = -VSmax * max([scr, 0]) / (KmS + scr);
lb(2,4) = -VmaxS * max([scr, 0]) / (KmS + scr + scr^2/20000) * (1 - eth/83);
ub(2,4) = Inf;

lb(2,5) = 0;
ub(2,5) = Inf;

lb(2,6) = -VmaxO * max([oxy, 0]) / (KmO + oxy);
ub(2,6) = Inf;

lb(2,7) = -VmaxPho * max([pho, 0]) / (KmPho + pho);
ub(2,7) = Inf;

lb(2,8) = -VmaxSul * max([sul, 0]) / (KmSul + sul);
ub(2,8) = Inf;

% Cheat - can uptake glucose and fructose but not sucrose
lb(3,1) = 0;
ub(3,1) = Inf;

lb(3,2) = -VmaxG * max([glc, 0]) / (KmG + glc);
ub(3,2) = Inf;

lb(3,3) = -VmaxF * max([fru, 0]) / (KmF + fru);
ub(3,3) = Inf;

lb(3,4) = 0;
ub(3,4) = Inf;

lb(3,5) = 0;
ub(3,5) = Inf;

lb(3,6) = -VmaxO * max([oxy, 0]) / (KmO + oxy);
ub(3,6) = Inf;

lb(3,7) = 0;
ub(3,7) = Inf;

lb(3,8) = 0;
ub(3,8) = Inf;


