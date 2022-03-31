
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

%% Assign value
N     = INFO.N;
ns    = INFO.ns;
param = INFO.param;

VGmax = param(1);
KmG   = param(2);
VFmax = param(3);
KmF   = param(4);
VSmax = param(5);
KmS   = param(6);
VOmax = param(7);
KmO   = param(8);

%% Set extracellular metabolite concentrations
glc = y(3);
fru = y(4);
scr = y(5);
oxy = y(8);

for j = 1:N
    % Private – can uptake glucose, fructose and sucrose
    lb(j,1) = 0;
    ub(j,1) = Inf;
    
    lb(j,2) = -VGmax * max([glc(j), 0]) / (KmG + glc(j));
    ub(j,2) = Inf;

    lb(j,3) = -VFmax * max([fru(j), 0]) / (KmF + fru(j));
    ub(j,3) = Inf;
    
    lb(j,4) = -VSmax * max([scr(j), 0]) / (KmS + scr(j));
    ub(j,4) = Inf;
        
    lb(j,5) = 0;
    ub(j,5) = Inf;
    
    lb(j,6) = -VOmax * max([oxy(j), 0]) / (KmO + oxy(j));
    ub(j,6) = Inf;
    
    % Cheat - can uptake glucose and fructose but not sucrose
    lb(N+j,1) = 0;
    ub(N+j,1) = Inf;
    
    lb(N+j,2) = -VGmax * max([glc(j), 0]) / (KmG + glc(j));
    ub(N+j,2) = Inf;

    lb(N+j,3) = -VFmax * max([fru(j), 0]) / (KmF + fru(j));
    ub(N+j,3) = Inf;
    
    lb(N+j,4) = 0;
    ub(N+j,4) = Inf;
        
    lb(N+j,5) = 0;
    ub(N+j,5) = Inf;
    
    lb(N+j,6) = -VOmax * max([oxy(j), 0]) / (KmO + oxy(j));
    ub(N+j,6) = Inf;
    
end



