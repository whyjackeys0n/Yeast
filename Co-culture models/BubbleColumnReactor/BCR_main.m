
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
clear
close
tic    
npoints = 21;
tspan = [0,1000];
nspecies = 2; 
nmodel = nspecies*npoints;

%% load GEM and manually set some bounds
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

% DFBAlab formality
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

%%
% Pass parameters
INFO.nmodel = nmodel;
INFO.DB = DB;
INFO.exID = exID;
INFO.C = C; 

N = npoints; 
ns = 11;
INFO.N = N;
INFO.ns = ns;

% Set operating conditions
L = 5;                   
L2D = 5/sqrt(3*4/pi);
Area = pi*((L/L2D)^2)/4;    
ug0 = 150/3600*(L/5);        
Gi = 200; 
ul = -50; 
db0 = 1.5;
zs = L/(N-1);                  
dl = .25;                     
tr = 310.15;                   
pL = 1.013e5;         
pc = 0.7;                  
Hc = 8.0e-4;    
ph = 0;          
Hh = 6.6e-4;  
pn = 0.3;
pc2 = 0;  
Hc2 = 2.5e-2;  
g = 9.81;     
densityL = 993.34;
viscosityL = 0.0009242;
surfaceTL = 40e-3;
D = 0.12;
Qmedia = D*L*Area;

% Calculate gas and liquid holdups
eg = 0.1;
el = 1-eg;
po = pL+1000*9.81*L*el;
cgi = pc*po/8.314/tr;           
hgi = ph*po/8.314/tr;       
c2gi = pc2*po/8.314/tr;       

% Set mass transfer coefficients
klc = 1e-4*3600;
klac = klc*6*eg/(1-eg)/db0*1000;
klah = klac;              
klac2 = klac;   

% Form condition vector
condit = [klac,Hc,pc,pc2,klac2,Hc2,tr,zs,N,ns,ug0,ul,dl,eg,el,cgi,pL,...
    c2gi,Area,Qmedia,D,densityL,g,viscosityL,db0,pn,klc,L,surfaceTL,Gi,ph,Hh];

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

% Pass parameters
INFO.param = param;
INFO.ns = ns;
INFO.N = N;
INFO.condit = condit;

% import steady-state solutions as initial conditions
load InitialConditions.mat
clsi = yfin3i;
cgi = yfin1i;
c2lsi = yfin4i;
c2gi = yfin2i;
% Set initial conditions
yo = [];
for i = 1:N
    r(i) = (cgi(i)+c2gi(i))/(cgi(1)+c2gi(1));
    eg(i) = 0.1;
    p(i) = pL + densityL*g*zs*(N-i)*(1-eg(i));
    rpc = (0.5-0.35)/(N-1)*(i-1);
    jg(i) = ug0*(pL/p(i))*r(i);
    db(i) = nthroot((pL/p(i))*r(i),3)*db0;
    ub(i) = 0.33*(g^0.76)*((densityL/viscosityL)^0.5)*((db(i)*1e-3).^1.28);
    yz = [cgi(1) 0 0 clsi(1) 0 0 p(i) jg(i) db(i) ub(i) eg(i)];
    yo = [yo yz];
end
yo = [yo, 0.01, 0.01, 0, 0, 0, 0, Gi, 0];
yo = [yo zeros(1,nmodel)];
Y0 = yo;

% Mass Matrix
M=[];
for i=1:N
    Mi = zeros(ns,ns);
    Mi(1,1) = 1;   
    Mi(2,2) = 1;  
    Mi(3,3) = 1;   
    Mi(4,4) = 1;   
    Mi(5,5) = 1;  
    Mi(6,6) = 1;   
    Mi(8,8) = 1; 
    M = blkdiag(M,Mi);
end

M = blkdiag(M,1);
M = blkdiag(M,1);
M = blkdiag(M,1);
M = blkdiag(M,1);
M = blkdiag(M,1);
M = blkdiag(M,1);
M = blkdiag(M,1);
M = blkdiag(M,1);

K = eye(nmodel,nmodel);
M = blkdiag(M,K);

% CPLEX Objects construction parameters
INFO.LPsolver = 1; % CPLEX = 0, Gurobi = 1.
                   % CPLEX works equally fine with both methods.
                   % Gurobi seems to work better with Method = 1, and 
                   % Mosek with Method = 0.
INFO.tol = 1E-7; % Feasibility, optimality and convergence tolerance for Cplex (tol>=1E-9). 
                 % It is recommended it is at least 2 orders of magnitude
                 % tighter than the integrator tolerance. 
                 % If problems with infeasibility messages, tighten this
                 % tolerance.
INFO.tolPh1 = INFO.tol; % Tolerance to determine if a solution to phaseI equals zero.
                   % It is recommended to be the same as INFO.tol. 
INFO.tolevt = 10*INFO.tol;%2E-7; %10*INFO.tol; % Tolerance for event detection. Has to be greater 
                   % than INFO.tol.

% You can modify the integration tolerances here.
% If some of the flows become negative after running the simulation once
% you can add the 'Nonnegative' option.

% options = odeset('AbsTol',1E-3,'RelTol',1E-2, ...
% 'NonNegative',1:N*ns,'Events',@evts, ...
% 'OutputSel',1,'OutputFcn',@odeplot, ...
% 'JPattern',syngas_bubblecolumn5_jac(N-1,ns));

% options = odeset('AbsTol',1E-6,'RelTol',1E-6, ...
% 'NonNegative',1:(N*ns+8),'Events',@evts, ...
% 'OutputSel',1,'OutputFcn',@odeplot);
NN = 1:length(Y0);
options = odeset('AbsTol',1E-5,'RelTol',1E-5,'Events',@evts,'Mass',M);

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
    
% Look at MATLAB documentation if you want to change solver.
% ode15s is more or less accurate for stiff problems. 
%     [T,Y] = ode15s(@DRHS,tspan,Y0,options,INFO);
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


elapsedTime = toc;
disp(['End Time = ',num2str(elapsedTime)]);

T = TF;
Y = YF;
%%
coeff=[];
for k=1:length(T)
    [~, ~, coeff(k,:),ca_mu(k,:),er_mu(k,:),ca_c(k,:),ca_c2(k,:),ca_a(k,:),...
        ca_e(k,:),ca_bdo(k,:),ca_lac(k,:),ca_h2(k,:),er_g(k,:),er_c2(k,:)...
        ,er_but(k,:),er_a(k,:),er_lac(k,:),er_h2(k,:)]=DRHS(T(k),Y(k,:),INFO);
end   
%% Plotting

%----------------------The first node point--------------------------------
figure(1)
% 
subplot(3,4,1)
hold on
plot(T,Y(:,ns*N+1))
xlabel('Time [h]')
ylabel('Ca Biomass [g/L]')

subplot(3,4,2)
hold on
plot(T,Y(:,ns*N+2))
xlabel('Time [h]')
ylabel('Er Biomass [g/L]')

subplot(3,4,3)
hold on
plot(T,Y(:,4))
xlabel('Time [h]')
ylabel('CO in liquid [mmol/L]')

subplot(3,4,4)
hold on
plot(T,Y(:,1))
xlabel('Time [h]')
ylabel('CO in gas [mmol/L]')

subplot(3,4,5)
hold on
plot(T,Y(:,ns*N+3))
xlabel('Time [h]')
ylabel('Acetate [mmol/L]')

subplot(3,4,6)
hold on
plot(T,Y(:,ns*N+4))
xlabel('Time [h]')
ylabel('Ethanol [mmol/L]')

subplot(3,4,7)
hold on
plot(T,Y(:,ns*N+7))
xlabel('Time [h]')
ylabel('Glucose [mmol/L]')

subplot(3,4,8)
hold on
plot(T,Y(:,ns*N+8))
xlabel('Time [h]')
ylabel('Butyrate [mmol/L]')

subplot(3,4,9)
hold on
plot(T,Y(:,7))
xlabel('Time [h]')
ylabel('Pressure [Pa]')

subplot(3,4,10)
hold on
plot(T,Y(:,8)*3600)
xlabel('Time [h]')
ylabel('Superficial gas velocity [m/h]')

subplot(3,4,11)
hold on
plot(T,Y(:,9))
xlabel('Time [h]')
ylabel('Bubble diameter [mm]')

subplot(3,4,12)
hold on
plot(T,Y(:,11))
xlabel('Time [h]')
ylabel('Gas holdup [-]')

%--------------The middle node point---------------------------------------
mz = round(N/2);

figure(2)

subplot(3,4,1)
hold on
plot(T,Y(:,ns*N+1))
xlabel('Time [h]')
ylabel('Ca Biomass [g/L]')

subplot(3,4,2)
hold on
plot(T,Y(:,ns*N+2))
xlabel('Time [h]')
ylabel('Er Biomass [g/L]')

subplot(3,4,3)
hold on
plot(T,Y(:,ns*(mz-1)+4))
xlabel('Time [h]')
ylabel('CO in liquid [mmol/L]')

subplot(3,4,4)
hold on
plot(T,Y(:,ns*(mz-1)+1))
xlabel('Time [h]')
ylabel('CO in gas [mmol/L]')

subplot(3,4,5)
hold on
plot(T,Y(:,ns*N+3))
xlabel('Time [h]')
ylabel('Acetate [mmol/L]')

subplot(3,4,6)
hold on
plot(T,Y(:,ns*N+4))
xlabel('Time [h]')
ylabel('Ethanol [mmol/L]')

subplot(3,4,7)
hold on
plot(T,Y(:,ns*N+7))
xlabel('Time [h]')
ylabel('Glucose [mmol/L]')

subplot(3,4,8)
hold on
plot(T,Y(:,ns*N+8))
xlabel('Time [h]')
ylabel('Butyrate [mmol/L]')

subplot(3,4,9)
hold on
plot(T,Y(:,ns*(mz-1)+7))
xlabel('Time [h]')
ylabel('Pressure [Pa]')

subplot(3,4,10)
hold on
plot(T,Y(:,ns*(mz-1)+8))
xlabel('Time [h]')
ylabel('Superficial gas velocity [m/s]')

subplot(3,4,11)
hold on
plot(T,Y(:,ns*(mz-1)+9))
xlabel('Time [h]')
ylabel('Bubble diameter [mm]')

subplot(3,4,12)
hold on
plot(T,Y(:,ns*(mz-1)+11))
xlabel('Time [h]')
ylabel('Gas holdup [-]')

%-------------------Outlet-------------------------------------------------
figure(3)

subplot(3,4,1)
hold on
plot(T,Y(:,ns*N+1))
xlabel('Time [h]')
ylabel('Ca Biomass [g/L]')

subplot(3,4,2)
hold on
plot(T,Y(:,ns*N+2))
xlabel('Time [h]')
ylabel('Er Biomass [g/L]')

subplot(3,4,3)
hold on
plot(T,Y(:,ns*(N-1)+4))
xlabel('Time [h]')
ylabel('CO in liquid [mmol/L]')

subplot(3,4,4)
hold on
plot(T,Y(:,ns*(N-1)+1))
xlabel('Time [h]')
ylabel('CO in gas [mmol/L]')

subplot(3,4,5)
hold on
plot(T,Y(:,ns*N+3))
xlabel('Time [h]')
ylabel('Acetate [mmol/L]')

subplot(3,4,6)
hold on
plot(T,Y(:,ns*N+4))
xlabel('Time [h]')
ylabel('Ethanol [mmol/L]')

subplot(3,4,7)
hold on
plot(T,Y(:,ns*N+7))
xlabel('Time [h]')
ylabel('Glucose [mmol/L]')

subplot(3,4,8)
hold on
plot(T,Y(:,ns*N+8))
xlabel('Time [h]')
ylabel('Butyrate [mmol/L]')

subplot(3,4,9)
hold on
plot(T,Y(:,ns*(N-1)+7))
xlabel('Time [h]')
ylabel('Pressure [Pa]')

subplot(3,4,10)
hold on
plot(T,Y(:,ns*(N-1)+8)*3600)
xlabel('Time [h]')
ylabel('Superficial gas velocity [m/h]')

subplot(3,4,11)
hold on
plot(T,Y(:,ns*(N-1)+9))
xlabel('Time [h]')
ylabel('Bubble diameter [mm]')

subplot(3,4,12)
hold on
plot(T,Y(:,ns*(N-1)+11))
xlabel('Time [h]')
ylabel('Gas holdup [-]')

%---------------spatial profiles versus time-------------------------------
for i=1:(N-1)
    zt(i)=i*zs;
end

klac = zeros(1, N);

for i=1:N
    yfin1(i)=Y(end,ns*(i-1)+1);
    yfin2(i)=Y(end,ns*(i-1)+2);
    yfin3(i)=Y(end,ns*(i-1)+3);
    yfin4(i)=Y(end,ns*(i-1)+4);
    yfin5(i)=Y(end,ns*(i-1)+5);
    yfin6(i)=Y(end,ns*(i-1)+6);
    yfin7(i)=Y(end,ns*(i-1)+7);
    yfin8(i)=Y(end,ns*(i-1)+8);
    yfin9(i)=Y(end,ns*(i-1)+9);
    yfin10(i)=Y(end,ns*(i-1)+10);
    yfin11(i)=Y(end,ns*(i-1)+11);
    klac(i) = 6*Y(end,ns*(i-1)+11)./(1-Y(end,ns*(i-1)+11))./(Y(end,ns*(i-1)+9)/1000)*klc;
    yfincabiomass(i) = Y(end,ns*N+1);
    yfinerbiomass(i) = Y(end,ns*N+2);
    yfinA(i) = Y(end,ns*N+3);
    yfinE(i) = Y(end,ns*N+4);
    yfinB(i) = Y(end,ns*N+5);
    yfinL(i) = Y(end,ns*N+6);
    yfinG(i) = Y(end,ns*N+7);
    yfinBUT(i) = Y(end,ns*N+8);
end

zt = [0 zt];

figure(4)
subplot(3,4,1)
hold on
plot(zt,yfincabiomass)
xlabel('Location [m]')
ylabel('Ca Biomass [g/L]')
xlim([0 L])

subplot(3,4,2)
hold on
plot(zt,yfinerbiomass)
xlabel('Location [m]')
ylabel('Er Biomass [g/L]')
xlim([0 L])

subplot(3,4,3)
hold on
plot(zt,yfin4)
xlabel('Location [m]')
ylabel('CO in liquid [mmol/L]')
xlim([0 L])

subplot(3,4,4)
hold on
plot(zt,yfin1)
xlabel('Location [m]')
ylabel('CO in gas [mmol/L]')
xlim([0 L])

subplot(3,4,5)
hold on
plot(zt,yfinA)
xlabel('Location [m]')
ylabel('Acetate [mmol/L]')
xlim([0 L])

subplot(3,4,6)
hold on
plot(zt,yfinE)
xlabel('Location [m]')
ylabel('Ethanol [mmol/L]')
xlim([0 L])

subplot(3,4,7)
hold on
plot(zt,yfinG)
xlabel('Location [m]')
ylabel('Glucose [mmol/L]')
xlim([0 L])

subplot(3,4,8)
hold on
plot(zt,yfinBUT)
xlabel('Location [m]')
ylabel('Butyrate [mmol/L]')
xlim([0 L])

subplot(3,4,9)
hold on
plot(zt,yfin7)
xlabel('Location [m]')
ylabel('Pressure [Pa]')
xlim([0 L])

subplot(3,4,10)
hold on
plot(zt,yfin8)
xlabel('Location [m]')
ylabel('Superficial gas velocity [m/s]')
xlim([0 L])

subplot(3,4,11)
hold on
plot(zt,yfin9)
xlabel('Location [m]')
ylabel('Bubble diameter [mm]')
xlim([0 L])

subplot(3,4,12)
hold on
plot(zt,yfin11)
xlabel('Location [m]')
ylabel('Gas holdup [-]')
xlim([0 L])

%---------------spatial profiles versus time (flow rate)-------------------
% Set other constants
Ma = 60.0/1000;   
Me = 46/1000; 
Mbdo = 90/1000;   
Ml = 90/1000; 
Mg = 180/1000;     
Mbut = 87/1000;   
% transfer concentration into mass flow rate
for i=1:N
     Y(:,ns*(i-1)+1) = Y(:,ns*(i-1)+1).*Area.*Y(:,ns*(i-1)+8);
     Y(:,ns*(i-1)+2) = Y(:,ns*(i-1)+2).*Area.*Y(:,ns*(i-1)+8);
     Y(:,ns*(i-1)+3) = Y(:,ns*(i-1)+3).*Area.*Y(:,ns*(i-1)+8);
end
vl_flowrate = -ul*Area/3600;
for i=1:N
     Y(:,ns*(i-1)+4) = Y(:,ns*(i-1)+4)*vl_flowrate;
     Y(:,ns*(i-1)+5) = Y(:,ns*(i-1)+5)*vl_flowrate;
     Y(:,ns*(i-1)+6) = Y(:,ns*(i-1)+6)*vl_flowrate;
end

yfincabiomass = yfincabiomass*L*1000*Area*D/1000;
yfinerbiomass = yfinerbiomass*L*1000*Area*D/1000; 
yfinA = yfinA*Ma*L*1000*Area*D/1000;
yfinE = yfinE*Me*L*1000*Area*D/1000;
yfinB = yfinB*Mbdo*L*1000*Area*D/1000;
yfinL = yfinL*Ml*L*1000*Area*D/1000;
yfinG = yfinG*Mg*L*1000*Area*D/1000;
yfinBUT = yfinBUT*Mbut*L*1000*Area*D/1000;

for i=1:(N-1)
    zt_2(i)=i*zs;
end

for i=1:N
    yfin1_2(i)=Y(end,ns*(i-1)+1);
    yfin2_2(i)=Y(end,ns*(i-1)+2);
    yfin3_2(i)=Y(end,ns*(i-1)+3);
    yfin4_2(i)=Y(end,ns*(i-1)+4);
    yfin5_2(i)=Y(end,ns*(i-1)+5);
    yfin6_2(i)=Y(end,ns*(i-1)+6);
    yfin7_2(i)=Y(end,ns*(i-1)+7);
    yfin8_2(i)=Y(end,ns*(i-1)+8);
    yfin9_2(i)=Y(end,ns*(i-1)+9);
    yfin10_2(i)=Y(end,ns*(i-1)+10);
    yfin11_2(i)=Y(end,ns*(i-1)+11);
end

zt_2 = [0 zt_2];
figure(5)
subplot(3,4,1)
hold on
plot(zt_2,yfincabiomass)
xlabel('Location [m]')
ylabel('Ca Biomass [kg/h]')
xlim([0 L])

subplot(3,4,2)
hold on
plot(zt_2,yfinerbiomass)
xlabel('Location [m]')
ylabel('Er Biomass [kg/h]')
xlim([0 L])

subplot(3,4,3)
hold on
plot(zt_2,yfin4_2)
xlabel('Location [m]')
ylabel('CO in liquid [mol/s]')
xlim([0 L])

subplot(3,4,4)
hold on
plot(zt_2,yfin1_2)
xlabel('Location [m]')
ylabel('CO in gas [mol/s]')
xlim([0 L])

subplot(3,4,5)
hold on
plot(zt_2,yfinA)
xlabel('Location [m]')
ylabel('Acetate [kg/h]')
xlim([0 L])

subplot(3,4,6)
hold on
plot(zt_2,yfinE)
xlabel('Location [m]')
ylabel('Ethanol [kg/h]')
xlim([0 L])

subplot(3,4,7)
hold on
plot(zt,yfinG)
xlabel('Location [m]')
ylabel('Glucose [kg/h]')
xlim([0 L])

subplot(3,4,8)
hold on
plot(zt,yfinBUT)
xlabel('Location [m]')
ylabel('Butyrate [kg/h]')
xlim([0 L])

subplot(3,4,9)
hold on
plot(zt,yfin7_2)
xlabel('Location [m]')
ylabel('Pressure [Pa]')
xlim([0 L])

subplot(3,4,10)
hold on
plot(zt,yfin8_2)
xlabel('Location [m]')
ylabel('Superficial gas velocity [m/s]')
xlim([0 L])

subplot(3,4,11)
hold on
plot(zt,yfin9_2)
xlabel('Location [m]')
ylabel('Bubble diameter [mm]')
xlim([0 L])

subplot(3,4,12)
hold on
plot(zt,yfin11_2)
xlabel('Location [m]')
ylabel('Gas holdup [-]')
xlim([0 L])

figure(6)
for i=1:N
    ca_mu_spatial(i)=ca_mu(end,i);
    er_mu_spatial(i)=er_mu(end,i);
    ca_c_spatial(i)=ca_c(end,i);
    ca_a_spatial(i)=ca_a(end,i);
    ca_e_spatial(i)=ca_e(end,i);
    ca_h_spatial(i) = ca_h2(end,i);
    er_g_spatial(i)=er_g(end,i);
    er_but_spatial(i)=er_but(end,i);
    er_a_spatial(i)=er_a(end,i);
    er_h_spatial(i)=er_h2(end,i);
end

subplot(3,3,1)
hold on
plot(zt,ca_mu_spatial)
xlabel('Location [m]')
ylabel('Ca Growth rate [1/h]')
xlim([0 L])

subplot(3,3,2)
hold on
plot(zt,er_mu_spatial)
xlabel('Location [m]')
ylabel('Er Growth rate [1/h]')
xlim([0 L])

subplot(3,3,3)
hold on
plot(zt,ca_c_spatial)
xlabel('Location [m]')
ylabel('Ca CO uptake rate [mmol/g/h]')
xlim([0 L])

subplot(3,3,4)
hold on
plot(zt,ca_a_spatial)
xlabel('Location [m]')
ylabel('Ca acetate synthesis rate [mmol/g/h]')
xlim([0 L])

subplot(3,3,5)
hold on
plot(zt,ca_e_spatial)
xlabel('Location [m]')
ylabel('Ca ethanol synthesis rate [mmol/g/h]')
xlim([0 L])

subplot(3,3,6)
hold on
plot(zt,er_g_spatial)
xlabel('Location [m]')
ylabel('Er glucose uptake rate [mmol/g/h]')
xlim([0 L])

subplot(3,3,7)
hold on
plot(zt,er_but_spatial)
xlabel('Location [m]')
ylabel('Er butyrate synthesis rate [mmol/g/h]')
xlim([0 L])

subplot(3,3,8)
hold on
plot(zt,er_a_spatial)
xlabel('Location [m]')
ylabel('Er acetate uptake rate [mmol/g/h]')
xlim([0 L])

subplot(3,3,9)
hold on
plot(zt,er_h_spatial)
xlabel('Location [m]')
ylabel('Er H2 synthesis rate [mmol/g/h]')
xlim([0 L])
