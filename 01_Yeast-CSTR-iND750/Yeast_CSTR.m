%% Script to simulate S. cerevisiae growth in a CSTR
% J. Wang, 1/14/2021

clear;

%% Set simulation time span and number of models
tic;
tspan = [0,10];    % unit: hour
nspecies = 1; 
npoints = 1;        % well-mixed : 1
nmodel = nspecies*npoints;

%% Load and setup extracellular variables
load yeastGEM 
modelwt = model;    % modelwt: model_wild_type
clear model

model = cell(1, npoints);   % pre-allocate the space for the loop
DB = zeros(1, npoints);     % pre-allocate the space for the loop
exID = cell(1, npoints);    % pre-allocate the space for the loop

for i = 1:npoints
    model{i} = modelwt; 
    DB(i) = 1000;
    iGr = find(strcmp(modelwt.rxnNames,'growth'));
    iGlu = find(strcmp(modelwt.rxnNames, 'D-glucose exchange'));
    iFru = find(strcmp(modelwt.rxnNames, 'D-fructose exchange'));
    iScr = find(strcmp(modelwt.rxnNames, 'sucrose exchange'));
    iO2 = find(strcmp(modelwt.rxnNames, 'oxygen exchange'));
    iEth = find(strcmp(modelwt.rxnNames, 'ethanol exchange'));
    exID{i} = [iGr iGlu iFru iScr iEth]; 
end

%% Setup lexicographic optimization objectives
minim = 1;  % defination for maximize and minimize
maxim = -1;

C = cell(1, npoints);   % pre-allocate the space for the loop

for i=1:npoints

    C{i}(1).sense = maxim;  % 1st objective
    C{i}(1).rxns = iGr;
    C{i}(1).wts = 1;

    C{i}(2).sense = minim;  % 2nd objective
    C{i}(2).rxns = iGlu;
    C{i}(2).wts = 1;

    C{i}(3).sense = minim;  % 3rd objective
    C{i}(3).rxns = iFru;
    C{i}(3).wts = 1;

    C{i}(4).sense = minim;  % 4th objective
    C{i}(4).rxns = iScr;
    C{i}(4).wts = 1;
    
    C{i}(5).sense = minim;  % 6th objective
    C{i}(5).rxns = iEth;
    C{i}(5).wts = 1;
    
end

%% Pass information to DFBAlab in INFO structure
INFO.nmodel = nmodel;   % nmodel = 1
INFO.DB = DB;           % DB = 1000
INFO.exID = exID;
INFO.C = C; 

%% Specify number of extracellular variables and pass information to DFBAlab
N = npoints;            % N = 1 
ns = length(exID{1});   % ns = 6
INFO.N = N;
INFO.ns = ns;

%% CSTR operation parameters and conditions
D = 0.1;    % dilution rate
Gf = 20;
Ff = 0;
Sf = 0;

%% Set uptake parameters
VGmax = 10;   % glucose
KmG = 0.5;  
VFmax = 10;   % fructose
KmF = 0.5;
VSmax = 0.5;  % sucrose
KmS = 0.5;
param = [VGmax; KmG; VFmax; KmF; VSmax; KmS];

%% Set initial conditions
Xi = 0.1;   % g/L
Gi = 20;    % mmol/L
Fi = 8;     % mmol/L
Si = 0.2;     % mmol/L
Ei = 0;     % mmol/L
yz = [Xi, Gi, Fi, Si, Ei];
Y0 = [yz, zeros(1,nmodel)];

%% Pass parameters
INFO.D = D;
INFO.Gf = Gf;
INFO.Ff = Ff;
INFO.Sf = Sf;
INFO.param = param;

%% Gurobi Objects construction parameters
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
options = odeset('AbsTol',1E-6,'RelTol',1E-6,'NonNegative',NN,'Events',@evts);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[model,INFO] = ModelSetupM(model,Y0,INFO);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

toc;

tic;

INFO.INFtimes = zeros(1,nmodel);
if INFO.LPsolver == 0
    [INFO] = LexicographicOpt(model,INFO);
elseif INFO.LPsolver == 1
    [INFO] = LexicographicOptG(model,INFO);
else
    disp('Solver not currently supported.');
end
toc;

tic;
tint = 0;
TF = [];
YF = [];
flux = [];
while tint<tspan(2)

    INFO.INFtimes;
    [T,Y] = ode15s(@DRHS_Yeast,tspan,Y0,options,INFO); 
    for i=1:length(T)
        [dy,fluxy] = DRHS_Yeast(T(i),Y(i,:),INFO);
        flux = [flux;fluxy];
    end
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
% Determine model with basis change
    value = evts(tint,Y0,INFO);
    ind = find(value(1:end-nmodel)<=0);
    fprintf('Basis change at time %d. ',tint);
    k = 0;
    ct = 0;
    
% Detect model infesibilities
    INFvector =  value(end-nmodel+1:end);
    for i=1:length(INFvector)
        if INFvector(i)>0 && INFO.INFtimes(i) == 0
            INFO.INFtimes(i) = tint;
            fprintf('Model %i is now dying at time %d. \n',i,tint);
        end
    end
    
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
               disp('Solver not currently supported.');
           end
           ind(ind2)=[];
        end
    end
end
toc;

T = TF;
Y = YF;

h=plot(T,Y(:,1),'-b',T,Y(:,2), '--r', T,Y(:,3),'-.',T,Y(:,4),'-g',T,Y(:,5));
set(h(1),'linewidth',2.5);
set(h(2),'linewidth',2.5);
set(h(3),'linewidth',2.5);
set(h(4),'linewidth',2.5);
set(h(5),'linewidth',2.5);
set(h(3),'Color',[0 0 0]);
legend('Biomass', 'Glucose', 'Fructose', 'Sucrose', 'Ethanol');
set(gca,'fontsize',24);
ylabel('Concentration [g/L]');
xlabel('Time [h]','FontSize',24)
ylim([0,25]);
