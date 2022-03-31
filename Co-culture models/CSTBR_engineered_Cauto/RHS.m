
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DFBAlab: Dynamic Flux Balance Analysis laboratory                       %
% Process Systems Engineering Laboratory, Cambridge, MA, USA              %
% July 2014                                                               %
% Written by Jose A. Gomez and Kai Höffner                                %
%                                                                         % 
% This code can only be used for academic purposes. When using this code  %
% please cite:                                                            %
%                                                                         %
% Gomez, J.A., Höffner, K. and Barton, P. I.                              %
% DFBAlab: A fast and reliable MATLAB code for Dynamic Flux Balance       %
% Analysis. Submitted.                                                    %
%                                                                         %
% COPYRIGHT (C) 2014 MASSACHUSETTS INSTITUTE OF TECHNOLOGY                %
%                                                                         %
% Read the LICENSE.txt file for more details.                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [lb,ub] = RHS(t,y,INFO)

% This subroutine updates the upper and lower bounds for the fluxes in the
% exID arrays in main. The output should be two matrices, lb and ub. The lb matrix
% contains the lower bounds for exID{i} in the ith row in the same order as
% exID. The same is true for the upper bounds in the ub matrix.
% Infinity can be used for unconstrained variables, however, it should be 
% fixed for all time. 

% Assign value

nmodel = INFO.nmodel;
ns = INFO.ns;
N = INFO.N;
param = INFO.param;

v_cm = param(1);
Kc   = param(2);
v_c2m = param(3);
Kc2   = param(4);
v_hm = param(5);
Kh   = param(6);
v_gm = param(7);
Kg = param(8);
v_am = param(9);
Ka = param(10);
Kic = param(11);
co_max = param(12);

j=1:N;
cl(j) = y(6);
c2l(j) = y(7);
hl(j) = y(8);
al(j) = y(9);
gl(j) = y(13);

for j = 1:N

    lb(j,1) = 0;
    ub(j,1) = Inf;
    
    lb(j,2) = -v_cm*max([cl(j) 0])/(Kc + cl(j) + cl(j)*cl(j)/Kic);
    ub(j,2) = 0;

    lb(j,3) = -v_c2m*max([c2l(j) 0])/(Kc2 + c2l(j));
    ub(j,3) = Inf;
    
    lb(j,4) = 0;
    ub(j,4) = Inf;
    
    lb(j,5) = 0;
    ub(j,5) = Inf;
    
    lb(j,6) = 0;
    ub(j,6) = Inf;
    
    lb(j,7) = 0;
    ub(j,7) = Inf;

    lb(j,8) = -v_hm*max([hl(j) 0])/(Kh + hl(j));
    ub(j,8) = 0;

    lb(j,9) = 0;
    ub(j,9) = Inf;

    lb(N+j,1) = 0;
    ub(N+j,1) = inf;

    lb(N+j,2) = (-v_gm*max([gl(j) 0])/(Kg + gl(j)))*(1-cl(j)/co_max);
    ub(N+j,2) = 0;

    lb(N+j,3) = 0;
    ub(N+j,3) = Inf;

    lb(N+j,4) = 0;
    ub(N+j,4) = Inf;

    lb(N+j,5) = -v_am*max([al(j) 0])/(Ka + al(j))*(1-cl(j)/co_max);
    ub(N+j,5) = 0;

    lb(N+j,6) = 0;
    ub(N+j,6) = Inf;

    lb(N+j,7) = 0;
    ub(N+j,7) = Inf;

end

end

