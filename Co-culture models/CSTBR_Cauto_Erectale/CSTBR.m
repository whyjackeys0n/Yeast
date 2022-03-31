clear
close
tic

tspan = [0,1000];
nspecies = 2;
npoints = 1;
nmodel = nspecies*npoints;

load Cauto_modelwt.mat     
load iEre400_norm.mat
modelwt.ub(960) = 1000;
model_er.lb(435:453) = 0;
model_er.lb(140:143) = -1000;
for i = 1:npoints
    model{i} = modelwt; 
    DB(i) = 1000;
    model{npoints+i} = model_er;
    DB(npoints+i) = 1000;
    exID{i}=[1026 960 854 972 969 986 658 966];
    exID{npoints+i} = [465 133 134 136 138 139 144];
end

minim = 1;
maxim = -1;

for i=1:npoints

    C{i}(1).sense = maxim;
    C{i}(1).rxns = [1026];
    C{i}(1).wts = [1];

    C{i}(2).sense = minim;
    C{i}(2).rxns = [960];
    C{i}(2).wts = [1];

    C{i}(3).sense = minim;
    C{i}(3).rxns = [854];
    C{i}(3).wts = [1];

    C{i}(4).sense = minim;
    C{i}(4).rxns = [972];
    C{i}(4).wts = [1];

    C{i}(5).sense = minim;
    C{i}(5).rxns = [969];
    C{i}(5).wts = [1];

    C{i}(6).sense = minim;
    C{i}(6).rxns = [986];
    C{i}(6).wts = [1];

    C{i}(7).sense = minim;
    C{i}(7).rxns = [658];
    C{i}(7).wts = [1];

    C{i}(8).sense = minim;
    C{i}(8).rxns = [966];
    C{i}(8).wts = [1];
    
    C{npoints+i}(1).sense = maxim;
    C{npoints+i}(1).rxns = [465];
    C{npoints+i}(1).wts = [1];

    C{npoints+i}(2).sense = minim;
    C{npoints+i}(2).rxns = [133];
    C{npoints+i}(2).wts = [1];

    C{npoints+i}(3).sense = minim;
    C{npoints+i}(3).rxns = [134];
    C{npoints+i}(3).wts = [1]; 

    C{npoints+i}(4).sense = minim;
    C{npoints+i}(4).rxns = [136];
    C{npoints+i}(4).wts = [1];

    C{npoints+i}(5).sense = minim;
    C{npoints+i}(5).rxns = [138];
    C{npoints+i}(5).wts = [1];

    C{npoints+i}(6).sense = minim;
    C{npoints+i}(6).rxns = [139];
    C{npoints+i}(6).wts = [1];

    C{npoints+i}(7).sense = minim;
    C{npoints+i}(7).rxns = [144];
    C{npoints+i}(7).wts = [1];
end

INFO.nmodel = nmodel;
INFO.DB = DB;
INFO.exID = exID;
INFO.C = C; 

N = npoints; 
ns = 14;
INFO.N = N;
INFO.ns = ns;

% operation parameters and conditions
tr = 310.15;                              
D = 0.12;                      
vg = 1;                    
V = 1;                        
eg = 0.1;
pL = 1.013e5;                     
pc = 0.7;                                 
Hc = 8.0e-4;                     
pc2 = 0.3;                              
Hc2 = 2.5e-2;             
ph = 0;                                
Hh = 6.6e-4;         

% Set mass transfer coefficients
klac = 100;                     
klac2 = klac;             
klah = klac;                  
el = 1-eg;

% Set uptake parameters
vcm = 50;  
Kmc = 0.1;
Kic = 5; 
vgm = 10;
Kmg = 0.5;
vam = 5;
Kma = 0.5;
co_max = 0.8;
param = [
     vcm                
     Kmc          
     vcm            
     Kmc         
     vcm               
     Kmc               
     vgm               
     Kmg             
     vam          
     Kma               
     Kic       
     co_max       
     ];

% Set initial conditions
yo = [];
Xcai = 0.01;
Xeri = 0.01;
cgi = pc*pL/8.314/tr;         
c2gi = pc2*pL/8.314/tr;        
hgi = ph*pL/8.314/tr;      
cli = 0;
c2li = 0;
hli = 0;
ai = 0;
ei = 0;
bdoi = 0;
li = 0;
gi = 10;
buti = 0;
yz = [Xcai, Xeri, cgi, c2gi, hgi, cli, c2li, hli, ai, ei, bdoi, li, gi, buti];
yo = [yo yz];
yo = [yo zeros(1,nmodel)];
Y0 = yo;

% Form condition vector
condit = [klac,Hc,klac2,Hc2,klah,Hh,tr,eg,el,cgi,c2gi,hgi,gi,D,vg];

% Pass parameters
INFO.param = param;
INFO.ns = ns;
INFO.N = N;
INFO.condit = condit;

% Gurobi Objects construction parameters
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
INFO.tolevt = 10*INFO.tol; % Tolerance for event detection. Has to be greater 
                   % than INFO.tol.

% You can modify the integration tolerances here.
% If some of the flows become negative after running the simulation once
% you can add the 'Nonnegative' option.

NN = 1:length(Y0);
options = odeset('AbsTol',1E-7,'RelTol',1E-7,'NonNegative',NN,'Events',@evts);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[model,INFO] = ModelSetupM(model,Y0,INFO);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
INFO.INFtimes = zeros(1,nmodel);
if INFO.LPsolver == 0
    [INFO] = LexicographicOpt(model,INFO);
elseif INFO.LPsolver == 1
    [INFO] = LexicographicOptG(model,INFO);
else
    display('Solver not currently supported.');
end

tint = 0;
TF = [];
YF = [];
flux = [];
while tint<tspan(2)

    INFO.INFtimes;
    [T,Y] = ode15s(@DRHS,tspan,Y0,options,INFO); 
    for i=1:length(T)
        [dy,fluxy] = DRHS(T(i),Y(i,:),INFO);
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
%Determine model with basis change
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
               display('Solver not currently supported.');
           end
           ind(ind2)=[];
        end
    end
end
T = TF;
Y = YF;

toc

%% Plot
figure(1)
subplot(3,4,1)
plot(T,Y(:,1))
xlabel('Time [h]')
ylabel('C.auto Biomass [g/L]')

subplot(3,4,2)
plot(T,Y(:,2))
xlabel('Time [h]')
ylabel('E.rectale Biomass [g/L]')

subplot(3,4,3)
plot(T,Y(:,3))
xlabel('Time [h]')
ylabel('CO gas concentration [mmol/L]')

subplot(3,4,4)
plot(T,Y(:,4))
xlabel('Time [h]')
ylabel('CO2 gas concentration [mmol/L]')

subplot(3,4,5)
plot(T,Y(:,5))
xlabel('Time [h]')
ylabel('H2 gas concentration [mmol/L]')

subplot(3,4,6)
plot(T,Y(:,6))
xlabel('Time [h]')
ylabel('CO liquid concentration [mmol/L]')

subplot(3,4,7)
plot(T,Y(:,7))
xlabel('Time [h]')
ylabel('CO2 liquid concentration [mmol/L]')

subplot(3,4,8)
plot(T,Y(:,8))
xlabel('Time [h]')
ylabel('H2 liquid concentration [mmol/L]')

subplot(3,4,9)
plot(T,Y(:,9))
xlabel('Time [h]')
ylabel('acetate concentration [mmol/L]')

subplot(3,4,10)
plot(T,Y(:,10))
xlabel('Time [h]')
ylabel('ethanol concentration [mmol/L]')

subplot(3,4,11)
plot(T,Y(:,13))
xlabel('Time [h]')
ylabel('glucose concentration [mmol/L]')

subplot(3,4,12)
plot(T,Y(:,14))
xlabel('Time [h]')
ylabel('butyrate concentration [mmol/L]')

%%

% Set other constants
Ma = 60.0/1000;     % acetate molecular weight g/mmol
Me = 46/1000;       % ethanol molecular weight g/mmol
Mbdo = 90/1000;     % 23BDO molecular weight g/mmol
Ml = 90/1000;     % Lactate molecular weight g/mmol
Mg = 180/1000;      % glucose molecular weight g/mmol
Mbut = 87/1000;     % butyrate molecular weight g/mmol

ca_mu = flux(:,1);
ca_c = flux(:,2);
ca_c2 = flux(:,3);
ca_a = flux(:,4);
ca_e = flux(:,5);
ca_bdo = flux(:,6);
ca_lac = flux(:,7);
ca_h = flux(:,8);
er_mu = flux(:,9);
er_g = flux(:,10);
er_c2 = flux(:,11);
er_but = flux(:,12);
er_a = flux(:,13);
er_lac = flux(:,14);
er_h = flux(:,15);

figure(2)
subplot(3,4,1)
plot(T,ca_mu(:,1))
xlabel('Time [h]')
ylabel('C.auto growth rate [1/h]')

subplot(3,4,2)
plot(T,er_mu(:,1))
ytickformat('%.3f')
xlabel('Time [h]')
ylabel('E.rectale growth rate [1/h]')

subplot(3,4,3)
plot(T,ca_c(:,1))
xlabel('Time [h]')
ylabel('CO uptake rate [mmol/gDW/h]')

subplot(3,4,4)
plot(T,ca_c2(:,1))
xlabel('Time [h]')
ylabel('CO2 uptake rate [mmol/gDW/h]')

subplot(3,4,5)
plot(T,ca_a(:,1))
xlabel('Time [h]')
ylabel('acetate secretion rate [mmol/gDW/h]')

subplot(3,4,6)
plot(T,ca_e(:,1))
xlabel('Time [h]')
ylabel('ethanol secretion rate [mmol/gDW/h]')

subplot(3,4,7)
plot(T,ca_h(:,1))
xlabel('Time [h]')
ylabel('H2 uptake rate [mmol/gDW/h]')

subplot(3,4,8)
plot(T,er_g(:,1))
xlabel('Time [h]')
ylabel('ER: glucose uptake rate [mmol/gDW/h]')

subplot(3,4,9)
plot(T,er_c2(:,1))
xlabel('Time [h]')
ylabel('ER: CO2 secretion rate')

subplot(3,4,10)
plot(T,er_but(:,1))
ytickformat('%.3f')
xlabel('Time [h]')
ylabel('ER: butyrate secretion rate')

subplot(3,4,11)
plot(T,er_a(:,1))
xlabel('Time [h]')
ylabel('ER: acetate uptake rate')

subplot(3,4,12)
plot(T,er_h(:,1))
xlabel('Time [h]')
ylabel('ER: H2 secretion rate')
