%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DFBAlab: Dynamic Flux Balance Analysis laboratory                       %
% Process Systems Engineering Laboratory, Cambridge, MA, USA              %
% July 2014                                                               %
% Written by Jose A. Gomez and Kai H�ffner                                %
%                                                                         % 
% This code can only be used for academic purposes. When using this code  %
% please cite:                                                            %
%                                                                         %
% Gomez, J.A., H�ffner, K. and Barton, P. I. (2014).                      %
% DFBAlab: A fast and reliable MATLAB code for Dynamic Flux Balance       %
% Analysis. BMC Bioinformatics, 15:409                                    % 
%                                                                         %
% COPYRIGHT (C) 2014 MASSACHUSETTS INSTITUTE OF TECHNOLOGY                %
%                                                                         %
% Read the LICENSE.txt file for more details.                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function dy = DRHS(t, y, INFO)

% Y1 = Volume (L)
% Y2 = Biomass EColi (gDW/L)
% Y3 = Glucose (g/L)
% Y4 = Xylose (g/L)
% Y5 = O2 (mmol/L)
% Y6 = Ethanol (g/L)
% Y7 = Penalty
Vol = y(1);
X(1) = y(2);
for i=1:4
   S(i) = y(2+i);
end

Fin = 0.1; %L/h
Fout = 0.1; %L/h

Xfeed(1) = 0;
Xfeed(2) = 0;

KhO2 = 1.3; % mmol/(L X atm)
k = 0.1; %1/h
Po2 = 0.21; %atm

MT(1) = 0;
MT(2) = 0; 
MT(3) = k*(Po2*KhO2-S(3));
MT(4) = 0;

MW(1) = 180.1559/1000;
MW(2) = 150.13/1000; 
MW(3) = 16/1000;
MW(4) = 46.06844/1000;

Sfeed(1) = 5; %g/L
Sfeed(2) = 2.5; %g/L
Sfeed(3) = Po2*KhO2; %mmol/L
Sfeed(4) = 0;


%% Update bounds and solve for fluxes
[flux,penalty] = solveModel(t,y,INFO);
%%

% The elements of the flux matrix have the sign given to them by the
% coefficients in the Cost vector in main. 
% Example, if:
% C{k}(i).rxns = [144, 832, 931];
% C{k}(i).wts = [3, 1, -1];
% Then the cost vector for this LP will be:
% flux(k,i) = 3*v_144 + v_832 - v_931 
% The penalty is an array containing the feasibility LP problem optimal
% objective function value for each model. 

v(1,4) = flux(1,2);
v(1,1) = flux(1,3);
v(1,2) = flux(1,4);
v(1,3) = flux(1,5);

%% Dynamics
dy = zeros(7,1);    % a column vector
dy(1) = Fin-Fout;    % Volume
dy(2) = flux(1,1)*X(1) + (Xfeed(1)*Fin - X(1)*Fout)/y(1);   % Biomass EColi
for i = 1:4
    dy(i+2) = v(1,i)*MW(i)*X(1) + MT(i) + (Sfeed(i)*Fin - S(i)*Fout)/y(1) ;  
end
dy(7) = penalty(1);
end