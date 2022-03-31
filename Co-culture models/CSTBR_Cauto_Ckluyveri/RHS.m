
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
v_am = param(7);
Ka = param(8);
v_em = param(9);
Ke = param(10);
v_sucm = param(11);
Ksuc = param(12);
v_prom = param(13);
Kpro = param(14);
v_crom = param(15);
Kcro = param(16);
v_vinylam = param(17);
Kvinyla = param(18);
v_propim = param(19);
Kpropi = param(20);
v_butm = param(21);
Kbut = param(22);
v_h2m = param(23);
Kh2 = param(24);
Kic = param(25);
co_max = param(26);

j=1:N;
cl(j) = y(6);
c2l(j) = y(7);
hl(j) = y(8);
al(j) = y(9);
el(j) = y(10);
sucl(j) = y(13);
prol(j) = y(14);
crol(j) = y(15);
vinylal(j) = y(16);
propil(j) = y(17);
butl(j) = y(18);

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

    lb(N+j,2) = -v_am*max([al(j) 0])/(Ka + al(j))*max([(1-cl(j)/co_max) 0]);
    ub(N+j,2) = Inf;

    lb(N+j,3) = -v_em*max([el(j) 0])/(Ke + el(j))*max([(1-cl(j)/co_max) 0]);
    ub(N+j,3) = Inf;

    lb(N+j,4) = -v_sucm*max([sucl(j) 0])/(Ksuc + sucl(j))*max([(1-cl(j)/co_max) 0]);
    ub(N+j,4) = Inf;

    lb(N+j,5) = -v_prom*max([prol(j) 0])/(Kpro + prol(j))*max([(1-cl(j)/co_max) 0]);
    ub(N+j,5) = 0;

    lb(N+j,6) = -v_crom*max([crol(j) 0])/(Kcro + crol(j))*max([(1-cl(j)/co_max) 0]);
    ub(N+j,6) = Inf;

    lb(N+j,7) = -v_vinylam*max([vinylal(j) 0])/(Kvinyla + vinylal(j))*max([(1-cl(j)/co_max) 0]);
    ub(N+j,7) = Inf;

    lb(N+j,8) = -v_propim*max([propil(j) 0])/(Kpropi + propil(j))*max([(1-cl(j)/co_max) 0]);
    ub(N+j,8) = Inf;

    lb(N+j,9) = -v_butm*max([butl(j) 0])/(Kbut + butl(j))*max([(1-cl(j)/co_max) 0]);
    ub(N+j,9) = Inf;

    lb(N+j,10) = 0;
    ub(N+j,10) = Inf;

    lb(N+j,11) = 0;
    ub(N+j,11) = Inf;

    lb(N+j,12) = -v_h2m*max([hl(j) 0])/(Kh2 + hl(j))*max([(1-cl(j)/co_max) 0]);
    ub(N+j,12) = Inf;    
end

end

