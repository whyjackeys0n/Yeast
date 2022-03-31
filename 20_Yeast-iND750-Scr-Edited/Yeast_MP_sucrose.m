clear;
clc;

%% Set simulation time span and number of models
tspan    = [0,24];  % unit: hour
nspecies = 3;       % public(1) & private(2) & cheat(3) 
npoints  = 1;       % well-mixed : 1
nmodel   = nspecies*npoints;

%% Load and setup extracellular variables
load iND750
load iND750_Private
load iND750_Cheat
modelpb = iND750;   % modelpb: model_wild_type
modelpr = iND750_Private;
modelch = iND750_Cheat;

iGr = find(strcmp(modelpb.rxnNames, 'biomass SC4 bal'));
iGlu = find(strcmp(modelpb.rxnNames, 'D-Glucose exchange'));
iFru = find(strcmp(modelpb.rxnNames, 'D-Fructose exchange'));
iScr = find(strcmp(modelpb.rxnNames, 'Sucrose exchange'));
iO2 = find(strcmp(modelpb.rxnNames, 'O2 exchange'));
iEth = find(strcmp(modelpb.rxnNames, 'Ethanol exchange'));
% iATP = find(strcmp(modelwt.rxnNames, 'ATP maintenance requirement'));
% modelwt.lb(iATP) = 10;
% modelwt.ub(iATP) = 1000;

model{1} = modelpb; 
DB(1) = 1000;
exID{1} = [iGr iGlu iFru iScr iEth iO2]; 

model{2} = modelpr; 
DB(2) = 1000;
exID{2} = [iGr iGlu iFru iScr iEth iO2]; 

model{3} = modelch; 
DB(3) = 1000;
exID{3} = [iGr iGlu iFru iScr iEth iO2]; 


%% Setup lexicographic optimization objectives
minim = 1;                  % defination for maximize 
maxim = -1;                 % defination for minimize

for i=1:nspecies
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
end

%% Pass information to DFBAlab in INFO structure
INFO.nmodel = nmodel;  
INFO.DB = DB;         
INFO.exID = exID;
INFO.C = C; 

%% Specify number of extracellular variables and pass to DFBAlab
ns = length(exID{1})*nmodel;
INFO.ns = ns;

%% CSTR operation parameters and conditions
D   = 0;        % dilution rate (h-1)
Dg  = 1000000;  % gas dilution rate (h-1)

P   = 101325;   % pressure (Pa) 
Tr  = 303.15;   % temperature (K)
R   = 8.314;    % gas constant (L·Pa·K-1·mmol-1)
Oc  = 0.21;     % oxygen mole fraction in feed gas
Kla = 100;      % volumetric gas-liquid mass transfer coefficient (h-1)

% Henry's Law coefficient (mol·L-1·atm-1)
% The Atmospheric Chemist’s Companion, 
% Numerical Data for Use in the Atmospheric Sciences, 
% Peter Warneck, Jonathan Williams, 2012, P290, Table 8.23
% ln(x) = A + B/T + ClnT
Oh  = exp(-175.33 + 8747.5/Tr + 24.453*log(Tr));

Gf  = 0;                % glucose feed (mmol·L-1)
Ff  = 0;                % fructose feed (mmol·L-1)
Sf  = 0;                % sucrose feed (mmol·L-1)
Ogf = Oc*P/R/Tr;        % oxygen feed in gas (mmol·L-1)

cond = [P, Tr, R, Oc, Kla, Oh];
feed = [Gf, Ff, Sf, Ogf];

%% Set uptake parameters
VGmax = 20;             % maximum glucose uptake rate (mmol/gDW/h)
KmG   = 2;              % glucose uptake saturation constants (mmol/L)
VFmax = VGmax;          % maximum fructose uptake rate (mmol/gDW/h)
KmF   = 2;              % fructose uptake saturation constants (mmol/L)
VSmax = 1;              % maximum sucrose uptake rate (mmol/gDW/h)
KmS   = 1;              % sucrose uptake saturation constants (mmol/L)
VOmax = 8;              % maximum oxygen uptake rate (mmol/gDW/h)
KmO   = 0.0031;         % oxygen uptake saturation constants (mmol/L)
param = [VGmax; KmG; VFmax; KmF; VSmax; KmS; VOmax; KmO];

%% Set enzyme synthesis and sucrose degradation model parameters
Alphah = 0.05;          % Basal enzyme synthesis rate
Betahm = 0.005;         % Maximum rate of induced synthesis (h-1)
Kbh    = 0.02;          % Saturation constant for induction (mmol/L)
Kibhh  = 0.07;          % Inhibition constant of glucose on induced beta glucosidase synthesis (mmol/L)
Gamah  = 0.3;           % Enzyme degradation rate constant (h-1)
Kp     = 0.1;           % Enzyme penalty constant

kdg    = 267.7;         % Reaction rate constant (h-1)
ks     = 25.5;          % Sucrose saturation constant (mmol/L)
kin    = 2.1;           % Inhibition constant for glucose (mmol/L)
enzyme = [Alphah; Betahm; Kbh; Kibhh; Gamah; Kp; kdg; ks; kin];

%% Pass conditions and parameters
INFO.D  = D;
INFO.Dg = Dg;
INFO.cond   = cond;
INFO.feed   = feed;
INFO.param  = param;
INFO.enzyme = enzyme;

%% Set initial conditions
Xpbi = 250*1e6*4.6e-11;     % g/L
Xpri = 250*1e6*4.6e-11;     % g/L
Xchi = 250*1e6*4.6e-11;     % g/L
Gi   = 0;                   % mmol/L
Fi   = 0;                   % mmol/L
Si   = 1*10/342*1000;  % mmol/L
Ei   = 0;                   % mmol/L
Ogi  = Ogf;                 % mmol/L <= 1*Pa/(L·Pa·K-1·mmol-1)/K
Oli  = Ogi*R*Tr/101325*Oh*1000;
% (mmol·L-1) <= (mmol·L-1)*(L·Pa·K-1·mmol-1)*K*(atm·Pa-1)*(mol·L-1·atm-1)*(mmol/mol)
Invi = 0;

yz = [Xpbi, Xpri, Xchi, Gi, Fi, Si, Ei, Ogi, Oli, Invi];
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
figure(1)
plot(T,Y(:,1),'linewidth',2.5);
hold on
plot(T,Y(:,2),'linewidth',2.5);
hold on
plot(T,Y(:,3),'linewidth',2.5);

legend('Public','Private','Cheat');
set(gca,'fontsize',24);
ylabel('Concentration [g(mmol)/L]');
xlabel('Time [h]','FontSize',24);
xlim([0,T(end)]);

% figure(2)
% plot(T,Y(:,4),'linewidth',2.5);
% hold on
% plot(T,Y(:,5),'linewidth',2.5);
% hold on
% plot(T,Y(:,6),'linewidth',2.5);
% 
% legend('Glucose','Fructose','Sucrose');
% set(gca,'fontsize',24);
% ylabel('Concentration [g(mmol)/L]');
% xlabel('Time [h]','FontSize',24);
% xlim([0,T(end)]);


