%% Script to simulate S. cerevisiae growth in a CSTR (Private & Cheat)
% J. Wang, 2/28/2021
clear;
clc;

%% Set simulation time span and number of models
tspan    = [0,50];  % unit: hour
nspecies = 2;       % private(1) & cheat(2) 
npoints  = 1;       % well-mixed : 1
nmodel   = nspecies*npoints;

%% Load and setup extracellular variables
load iND750
modelwt = iND750;   % modelwt: model_wild_type

iGr = find(strcmp(modelwt.rxnNames, 'biomass SC4 bal'));
iGlu = find(strcmp(modelwt.rxnNames, 'D-Glucose exchange'));
iFru = find(strcmp(modelwt.rxnNames, 'D-Fructose exchange'));
iScr = find(strcmp(modelwt.rxnNames, 'Sucrose exchange'));
iO2 = find(strcmp(modelwt.rxnNames, 'O2 exchange'));
iEth = find(strcmp(modelwt.rxnNames, 'Ethanol exchange'));

for i = 1:npoints
    model{i} = modelwt; 
    DB(i) = 1000;
    exID{i} = [iGr iGlu iFru iScr iEth iO2]; 
    
    model{npoints+i} = modelwt; 
    DB(npoints+i) = 1000;
    exID{npoints+i} = [iGr iGlu iFru iScr iEth iO2]; 
end

%% Setup lexicographic optimization objectives
minim = 1;  % defination for maximize and minimize
maxim = -1;

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
    
    C{i}(5).sense = minim;  % 5th objective
    C{i}(5).rxns = iEth;
    C{i}(5).wts = 1;
    
    C{i}(6).sense = minim;  % 6th objective
    C{i}(6).rxns = iO2;
    C{i}(6).wts = 1;
    
    C{npoints+i}(1).sense = maxim;  % 1st objective
    C{npoints+i}(1).rxns = iGr;
    C{npoints+i}(1).wts = 1;
    
    C{npoints+i}(2).sense = minim;  % 2nd objective
    C{npoints+i}(2).rxns = iGlu;
    C{npoints+i}(2).wts = 1;

    C{npoints+i}(3).sense = minim;  % 3rd objective
    C{npoints+i}(3).rxns = iFru;
    C{npoints+i}(3).wts = 1;

    C{npoints+i}(4).sense = minim;  % 4th objective
    C{npoints+i}(4).rxns = iScr;
    C{npoints+i}(4).wts = 1;
    
    C{npoints+i}(5).sense = minim;  % 5th objective
    C{npoints+i}(5).rxns = iEth;
    C{npoints+i}(5).wts = 1;
    
    C{npoints+i}(6).sense = minim;  % 6th objective
    C{npoints+i}(6).rxns = iO2;
    C{npoints+i}(6).wts = 1;
end

%% Pass information to DFBAlab in INFO structure
INFO.nmodel = nmodel;  
INFO.DB = DB;         
INFO.exID = exID;
INFO.C = C; 

%% Specify number of extracellular variables and pass information to DFBAlab
N = npoints; 
ns = length(exID{1})*nmodel;
INFO.N = N;
INFO.ns = ns;

%% CSTR operation parameters and conditions
D   = 0.25;     % dilution rate (h-1)
Dg  = 10;       % gas dilution rate (h-1)

P   = 101325;   % pressure (Pa) 
Tr  = 303.15;   % liquid temperature (K)
R   = 8.314;    % gas constant (L·Pa·K-1·mmol-1)
Oc  = 0.21;     % oxygen mole fraction in feed gas
Kla = 50;       % volumetric gas-liquid mass transfer coefficient (h-1)

% Henry's Law coefficient (mol·L-1·atm-1)
% The Atmospheric Chemist’s Companion, 
% Numerical Data for Use in the Atmospheric Sciences, 
% Peter Warneck, Jonathan Williams, 2012, P290, Table 8.23
% ln(x) = A + B/T + ClnT
Oh  = exp(-175.33 + 8747.5/Tr + 24.453*log(Tr));

Gf  = 20;       % glucose feed (mmol·L-1)
Ff  = 5;        % fructose feed (mmol·L-1)
Sf  = 1;        % sucrose feed (mmol·L-1)
Ogf = Oc*P/R/Tr;% oxygen feed in gas (mmol·L-1)

cond = [P, Tr, R, Oc, Kla, Oh];
feed = [Gf, Ff, Sf, Ogf];

%% Set uptake parameters
VGmax = 10;     % maximum glucose uptake rate (mmol/gDW/h)
KmG = 0.5;      % glucose uptake saturation constants (mmol/L)
VFmax = 10;     % maximum fructose uptake rate (mmol/gDW/h)
KmF = 0.5;      % fructose uptake saturation constants (mmol/L)
VSmax = 10;     % maximum sucrose uptake rate (mmol/gDW/h)
KmS = 0.5;      % sucrose uptake saturation constants (mmol/L)
VOmax = 8;      % maximum oxygen uptake rate (mmol/gDW/h)
KmO = 0.003;    % oxygen uptake saturation constants (mmol/L)
param = [VGmax; KmG; VFmax; KmF; VSmax; KmS; VOmax; KmO];

%% Pass conditions and parameters
INFO.D  = D;
INFO.Dg = Dg;
INFO.cond  = cond;
INFO.feed  = feed;
INFO.param = param;

%% Set initial conditions
Xpri = 1;       % g/L
Xchi = 1;       % g/L
Gi   = Gf;      % mmol/L
Fi   = Ff;      % mmol/L
Si   = Sf;      % mmol/L
Ei   = 0;       % mmol/L
Ogi  = Ogf;     % mmol/L <= 1*Pa/(L·Pa·K-1·mmol-1)/K
Oli  = Ogi*R*Tr/101325*Oh*1000;
% (mmol·L-1) <= (mmol·L-1)*(L·Pa·K-1·mmol-1)*K*(atm·Pa-1)*(mol·L-1·atm-1)*(mmol/mol)

yz = [Xpri, Xchi, Gi, Fi, Si, Ei, Ogi, Oli];
Y0 = [yz, zeros(1,nmodel)];

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

INFO.INFtimes = zeros(1,nmodel);
if INFO.LPsolver == 0
    [INFO] = LexicographicOpt(model,INFO);
elseif INFO.LPsolver == 1
    [INFO] = LexicographicOptG(model,INFO);
else
    disp('Solver not currently supported.');
end

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

T = TF;
Y = YF;

%% Plot
h = plot(T,Y(:,1),'-k',T,Y(:,2),'-',T,Y(:,3),'--b',T,Y(:,4),'--',T,Y(:,5),'--g',T,Y(:,6),':m',T,Y(:,7),'-.r',T,Y(:,8),'-.c');
set(h(1),'linewidth',2.5);
set(h(2),'linewidth',2.5,'color',[0.54118, 0.16863, 0.88627]);
set(h(3),'linewidth',2.5);
set(h(4),'linewidth',2.5,'color',[1, 0.64706, 0]);
set(h(5),'linewidth',2.5);
set(h(6),'linewidth',2.5);
set(h(7),'linewidth',2.5);
set(h(8),'linewidth',2.5);
legend('Biomass Private', 'Biomass Cheat', 'Glucose', 'Fructose', 'Sucrose', 'Ethanol','Oxygen (g)','Oxygen (l)');
set(gca,'fontsize',24);
ylabel('Concentration [g/L]');
xlabel('Time [h]','FontSize',24);
ylim([0,35]);
set(gcf,'position',[500,300,1500,1000]);
