%% Script to simulate S. cerevisiae growth in a CSTR (Private & Cheat)
% J. Wang, 3/15/2021
clear;
clc;

%% Set simulation time span and number of models
tspan    = [0,70];  % unit: hour
nspecies = 3;       % public(1) & private(2) & cheat(3) 
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
iATP = find(strcmp(modelwt.rxnNames, 'ATP maintenance requirement'));
% modelwt.lb(iATP) = 10;
% modelwt.ub(iATP) = 1000;

AAlb = -0.06;
% The “complete supplement mixture” contains the following amino acids (mg/L): 
% Adenine(10)-C5H5N5(135), L-Arginine(50)-C6H14N4O2(174), L-Aspartic acid(80)-C4H7NO4(133), 
iAdenine = find(strcmp(modelwt.rxnNames,'Adenine exchange'));
modelwt.lb(iAdenine) = AAlb;
iArginine = find(strcmp(modelwt.rxnNames,'L-Arginine exchange'));
modelwt.lb(iArginine) = AAlb;
iArginine = find(strcmp(modelwt.rxnNames,'L-Aspartate exchange'));
modelwt.lb(iArginine) = AAlb;

% L-Histidine HCl(20)-C6H9N3O2(155), L-Isoleucine(50)-C6H13NO2(131), L-Leucine(100)-C6H13NO2(131), 
iHistidine = find(strcmp(modelwt.rxnNames,'L-Histidine exchange'));
modelwt.lb(iHistidine) = AAlb;
iIsoleucine  = find(strcmp(modelwt.rxnNames,'L-Isoleucine exchange'));
modelwt.lb(iIsoleucine) = AAlb;
iLeucine = find(strcmp(modelwt.rxnNames,'L-Leucine exchange'));
modelwt.lb(iLeucine) = AAlb;

% L-Lysine HCl(50)-C6H14N2O2(146), L-Methionine(20)-C5H11O2NS(149), L-Phenylalanine(50)-C9H11NO2(165), 
iLysine = find(strcmp(modelwt.rxnNames,'L-Lysine exchange'));
modelwt.lb(iLysine) = AAlb;
iMethionine = find(strcmp(modelwt.rxnNames,'L-Methionine exchange'));
modelwt.lb(iMethionine) = AAlb;
iPhenylalanine = find(strcmp(modelwt.rxnNames,'L-Phenylalanine exchange'));
modelwt.lb(iPhenylalanine) = AAlb;

% L-Threonine(100)-C4H9NO3(119), L-Tryptophan(50)-C11H12N2O2(204), L-Tyrosine(50)-C9H11NO3(181), 
iThreonine = find(strcmp(modelwt.rxnNames,'L-Threonine exchange'));
modelwt.lb(iThreonine) = AAlb;
iTryptophan = find(strcmp(modelwt.rxnNames,'L-Tryptophan exchange'));
modelwt.lb(iTryptophan) = AAlb;
iTyrosine = find(strcmp(modelwt.rxnNames,'L-Tyrosine exchange'));
modelwt.lb(iTyrosine) = AAlb;

% Uracil(20)-C4H4N2O2(112), Valine(140)-C5H11NO2(117). 
iUracil = find(strcmp(modelwt.rxnNames,'Uracil exchange'));
modelwt.lb(iUracil) = AAlb;
iValine = find(strcmp(modelwt.rxnNames,'L-Valine exchange'));
modelwt.lb(iValine) = AAlb;

for i = 1:nspecies
    model{i} = modelwt; 
    DB(i) = 1000;
    exID{i} = [iGr iGlu iFru iScr iEth iO2]; 
end

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
Dg  = 1000;       % gas dilution rate (h-1)

P   = 101325;   % pressure (Pa) 
Tr  = 303.15;   % temperature (K)
R   = 8.314;    % gas constant (L·Pa·K-1·mmol-1)
Oc  = 0.21;     % oxygen mole fraction in feed gas
Kla = 18;      % volumetric gas-liquid mass transfer coefficient (h-1)

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
VGmax = 3.5/180*1000*0.36;   % maximum glucose uptake rate (mmol/gDW/h)
KmG = 0.5/180*1000;     % glucose uptake saturation constants (mmol/L)
VFmax = VGmax;          % maximum fructose uptake rate (mmol/gDW/h)
KmF = 4.6;              % fructose uptake saturation constants (mmol/L)
VSmax = VGmax/63.45*5.77;       % maximum sucrose uptake rate (mmol/gDW/h)
KmS = 5;                % sucrose uptake saturation constants (mmol/L)
VOmax = 8;              % maximum oxygen uptake rate (mmol/gDW/h)
KmO = 0.1/32;           % oxygen uptake saturation constants (mmol/L)
param = [VGmax; KmG; VFmax; KmF; VSmax; KmS; VOmax; KmO];

%% Pass conditions and parameters
INFO.D  = D;
INFO.Dg = Dg;
INFO.cond  = cond;
INFO.feed  = feed;
INFO.param = param;

%% Set initial conditions
Xpbi = 100*4.6e-5;      % g/L
Xpri = 0;               % g/L
Xchi = 0;               % g/L
Gi   = 10/180*1000;     % mmol/L
Fi   = 0;               % mmol/L
Si   = 0;               % mmol/L
Ei   = 0;               % mmol/L
Ogi  = Ogf;             % mmol/L <= 1*Pa/(L·Pa·K-1·mmol-1)/K
Oli  = Ogi*R*Tr/101325*Oh*1000;
% (mmol·L-1) <= (mmol·L-1)*(L·Pa·K-1·mmol-1)*K*(atm·Pa-1)*(mol·L-1·atm-1)*(mmol/mol)

yz = [Xpbi, Xpri, Xchi, Gi, Fi, Si, Ei, Ogi, Oli];
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
h1 = plot(T,Y(:,1),'-k');
set(h1(1),'linewidth',2.5);
legend('Biomass Public');
set(gca,'fontsize',24);
ylabel('Concentration [g(mmol)/L]');
xlabel('Time [h]','FontSize',24);
xlim([0,T(end)]);
hold on
% ylim([0,max(Y(:,1))]);

% figure(2)
% h1 = plot(T,Y(:,4),'r');
% set(h1(1),'linewidth',2.5);
% legend('Glucose');
% set(gca,'fontsize',24);
% ylabel('Concentration','FontSize',24);
% xlim([0,T(end)]);
% ylim([0,1.5*max([Y(:,1);Y(:,2);Y(:,3)])]);



% figure(4)
% h2 = plot(T,flux(:,2),T,flux(:,6));
% set(h2(1),'linewidth',2.5,'color',[0.54118, 0.16863, 0.88627]);
% set(h2(2),'linewidth',2.5,'color',[1, 0.64706, 0]);
% legend('Glucose', 'Oxygen(L)');
% set(gca,'fontsize',24);
% ylabel('Flux','FontSize',24);
% xlim([0,T(end)]);
% ylim([0,2]);

%% load data
load CENPK21C_0410_1pct_glucose_pbpr  
CENPK21C_0410_1pct_glucose_pbpr_pop_density = CENPK21C_0410_1pct_glucose_pbpr.CENPK21C_0410_1pct_glucose_pbpr_pop_density;
CENPK21C_0410_1pct_glucose_pbpr_time = CENPK21C_0410_1pct_glucose_pbpr.CENPK21C_0410_1pct_glucose_pbpr_time;
CENPK21C_0410_1pct_glucose_pbpr_conc_pct = CENPK21C_0410_1pct_glucose_pbpr.CENPK21C_0410_1pct_glucose_pbpr_conc_pct;

%% unit convert from cell/microliter to gram/liter
unit = 1e6 * 4.6e-11;

%% S288C_0410_1pct_glucose_pbprch
pop_density_gL = CENPK21C_0410_1pct_glucose_pbpr_pop_density.*unit;
time = CENPK21C_0410_1pct_glucose_pbpr_time;
conc_pct = CENPK21C_0410_1pct_glucose_pbpr_conc_pct;
species = {'Sim','Exp Public','Exp Private'};

%% figure 1
pop_density_gL_mean = [];
re_exp_num = 4;
species_num = 2;

for i = 1:species_num   
    pop_density_gL_mean_unit = zeros(1,length(pop_density_gL));
    for j = re_exp_num*(i-1)+1:re_exp_num*(i-1)+re_exp_num
        pop_density_gL_mean_unit = pop_density_gL_mean_unit + pop_density_gL(j,:);
    end
    pop_density_gL_mean = [pop_density_gL_mean;pop_density_gL_mean_unit./re_exp_num];
end
for i = 1:species_num
    plot(time,pop_density_gL_mean(i,:),'LineWidth',2); 
    %plot(time,log10(pop_density_gL_mean(i,:)),'LineWidth',2);    
    hold on
end
legend(species,'FontSize',24)
set(gca,'FontSize',24)

figure(2)
h1 = plot(T,Y(:,9),'b');
set(h1(1),'linewidth',2.5);
legend('Oxygen(L)');
set(gca,'fontsize',24);
ylabel('Concentration','FontSize',24);
xlim([0,T(end)]);

figure(3)
h1 = plot(T,Y(:,7),'r');
set(h1(1),'linewidth',2.5);
legend('Ethanol');
set(gca,'fontsize',24);
ylabel('Concentration','FontSize',24);
xlim([0,T(end)]);