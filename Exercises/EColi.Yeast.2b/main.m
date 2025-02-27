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

% This example is based on Hanly, T.J., Henson, M.A. Dynamic flux balance
% modeling of microbial co-cultures for efficient batch fermentation of
% glucose and xylose mixtures. Biotechnology and Bioengineering
% 108(2),376-385 (2011).
% The genome-scale metabolic model used was iJR904 from Reed,J.L, Vo, T.D.,
% Schilling, C.H., Palsson, B.O. An expanded genome-scale model of
% Escherichia coli K-12 (iJR904 GSM/GPR). Genome Biology 4(9), 54 (2003). 

clear all
INFO.nmodel = 2; % Number of models
% Load models. These should be .mat files generated by the COBRA toolbox. 
% When generating these files using the COBRA toolbox, a big number is used
% as infinity. This number should be fed to the DB vector (Default bound).

load iJR904.mat 
model{1} = iJR904; 
DB(1) = 2000;

load iND750.mat 
model{2} = iND750; 
DB(2) = 1000;
INFO.DB = DB;
%% exID array
% You can either search the reaction names by name or provide them directly
% in the exID array.
% RxnNames = {'EX_glc(e)', 'EX_ac(e)', 'biomass'};
% for i = 1:length(RxnNames)
%    [a,exID(i)] = ismember(RxnNames(i),model.rxns);
% end

exID{1}=[344, 429, 392, 329];
exID{2}=[428,458,407,420];
INFO.exID = exID;
% Lower bounds and upper bounds for these reactions should be provided in
% the RHS code. 

% This codes solves the LPs in standard form. Bounds on exchange fluxes in 
% the exID array can be modified directly on the first 2*n rows where n is 
% the number of exchange fluxes. Order will be lower bound, upper bound, 
% lower bound, upper bound in the same order as exID. 
%
% NOTE: All bounds on fluxes in the exID arrays are relaxed to -Inf and 
% + Inf. These bounds need to be updated if needed in the RHS file.

%% Cost vectors
% Usually the first cost vector will be biomass maximization, but it can
% be any other objective. The CPLEX objects will minimize by default. 
% Report only nonzero elements. 
% The structure should be:
% C{model} = struct
% Element C{k}(i) of C, is the cost structure i for model k. 
% C{k}(i).sense = +1 for minimize, or -1 for maximize.
% C{k}(i).rxns = array containing the reactions in this objective. 
% C{k}(i).wts = array containing coefficients for reactions reported in 
% rxns. Both arrays should have the same length. 
% Example, if:
% C{k}(i).rxns = [144, 832, 931];
% C{k}(i).wts = [3, 1, -1];
% Then the cost vector for this LP will be:
% Cost{k}(i) = 3*v_144 + v_832 - v_931 (fluxes for model k). 
% This cost vector will be either maximized or minimized depending on the
% value of C{k}(i).sense.

% In SBML files, usually production fluxes are positive and uptake fluxes
% are negative. Keep in mind that maximizing a negative flux implies 
% minimizing its absolute value.
% Different models can have different number of objectives. 

minim = 1;
maxim = -1;

% Maximize growth
C{1}(1).sense = maxim;
C{1}(1).rxns = [150];
C{1}(1).wts = [1];
% Maximize ethanol
C{1}(2).sense = maxim;
C{1}(2).rxns = [329];
C{1}(2).wts = [1];
% Maximize glucose
C{1}(3).sense = maxim;
C{1}(3).rxns = [344];
C{1}(3).wts = [1];
% Maximize xylose
C{1}(4).sense = maxim;
C{1}(4).rxns = [429];
C{1}(4).wts = [1];
% Maximize oxygen
C{1}(5).sense = maxim;
C{1}(5).rxns = [392];
C{1}(5).wts = [1];

% Yeast
% Maximize growth
C{2}(1).sense = maxim;
C{2}(1).rxns = [1266];
C{2}(1).wts = [1];
% Glucose
C{2}(2).sense = maxim;
C{2}(2).rxns = [428];
C{2}(2).wts = [1];
% O2
C{2}(3).sense = maxim;
C{2}(3).rxns = [458];
C{2}(3).wts = [1];
% Ethanol
C{2}(4).sense = maxim;
C{2}(4).rxns = [420];
C{2}(4).wts = [1];

INFO.C = C;
% Initial conditions
% Y1 = Volume (L)
% Y2 = Biomass EColi (gDW/L)
% Y3 = Biomass Yeast (gdW/L)
% Y4 = Glucose (g/L)
% Y5 = Xylose (g/L)
% Y6 = O2 (mmol/L)
% Y7 = Ethanol (g/L)
% Y8 = Penalty
Y0 = [1 0.03 0.03 15.5 8 0.24 0 0]';

% Time of simulation
tspan = [0,20];

% CPLEX Objects construction parameters
INFO.LPsolver = 1; % CPLEX = 0, Gurobi = 1.
                   % CPLEX works equally fine with both methods.
                   % Gurobi seems to work better with Method = 1, and 
                   % Mosek with Method = 0.
INFO.tol = 1E-9; % Feasibility, optimality and convergence tolerance for Cplex (tol>=1E-9). 
                 % It is recommended it is at least 2 orders of magnitude
                 % tighter than the integrator tolerance. 
                 % If problems with infeasibility messages, tighten this
                 % tolerance.
INFO.tolPh1 = INFO.tol; % Tolerance to determine if a solution to phaseI equals zero.
                   % It is recommended to be the same as INFO.tol. 
INFO.tolevt = 2*INFO.tol; % Tolerance for event detection. Has to be greater 
                   % than INFO.tol.

% You can modify the integration tolerances here.
% If some of the flows become negative after running the simulation once
% you can add the 'Nonnegative' option.
NN = 1:length(Y0);
options = odeset('AbsTol',1E-6,'RelTol',1E-6,'Nonnegative',NN,'Events',@evts);

% INFO: You can use the INFO struct to pass parameters. Don't use any of 
% the names already declared or: INFO.t (carries time information), 
% INFO.ncost, INFO.lexID, INFO.LlexID, INFO.lbct, INFO.ubct, INFO.sense,
% INFO.b, INFO.pair. 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[model,INFO] = ModelSetupM(model,Y0,INFO);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if INFO.LPsolver == 0
    [INFO] = LexicographicOpt(model,INFO);
elseif INFO.LPsolver == 1
    [INFO] = LexicographicOptG(model,INFO);
else
    display('Solver not currently supported.');
end

tic
tint = 0;
TF = [];
YF = [];
while tint<tspan(2)
% Look at MATLAB documentation if you want to change solver.
% ode15s is more or less accurate for stiff problems. 
    [T,Y] = ode15s(@DRHS,tspan,Y0,options,INFO);
    TF = [TF;T];
    YF = [YF;Y];
    tint = T(end);
    tspan = [tint,tspan(2)];
    Y0 = Y(end,:);
    if tint == tspan(2)
        break;
    end
    
% Update b vector 
[INFO] = bupdate(tint,Y0,INFO);
%Determine model with basis change
    value = evts(tint,Y0,INFO);
    ind = find(value<=0);
    fprintf('Basis change at time %d. ',tint);
    k = 0;
    ct = 0;
    while ~isempty(ind)
        k = k + 1;
        ct = ct + size(model{k}.A,1);
        ind2 = find(ind<=ct);
        if ~isempty(ind2)
           INFO.flagbasis = k; 
           fprintf('Model %i. \n',k);
           % Perform lexicographic optimization
           if INFO.LPsolver == 0
               [INFO] = LexicographicOpt(model,INFO);
           elseif INFO.LPsolver == 1
               [INFO] = LexicographicOptG(model,INFO);
           else
               display('Solver not currently supported.');
           end
           ind(ind2)=[];
        end
    end
end
display(toc);
T = TF;
Y = YF;
% 
% % %% Plotting
figure(1)
h=plot(T,Y(:,2),'-b',T,Y(:,3),T,Y(:,4), '--r', T,Y(:,5),'-.',T,Y(:,6),T,Y(:,7));
set(h(1),'linewidth',2.5);
set(h(2),'linewidth',2.5);
set(h(3),'linewidth',2.5);
set(h(4),'linewidth',2.5);
set(h(5),'linewidth',2.5);
set(h(6),'linewidth',2.5);
set(h(3),'Color',[0 0 0]);
legend('E. Coli', 'Yeast', 'Glucose', 'Xylose','O_2','Ethanol')
set(gca,'fontsize',24)
ylabel('Concentration [g/L]');
xlabel('Time [h]','FontSize',24)
ylim([0,16]);
% xlim([0,8.1]);

figure(2)
left= 0.1;
bottom3=0.1;
mid = 0.47;
width=0.8;
height=0.8;

ax1 = axes('Position', [left bottom3 width height]);
h=plot(ax1,T,Y(:,8));
set(h(1),'linewidth',2.5);
set(gca,'fontsize',24)
ylabel('Penalty State');
xlabel('Time [h]');

