function [dy,fluxy]= DRHS_Yeast(t, y, INFO)

%% Assign values
nmodel = INFO.nmodel;
N      = INFO.N;
D      = INFO.D;
Dg     = INFO.Dg;
cond   = INFO.cond;
feed   = INFO.feed;

P   = cond(1); 
Tr  = cond(2); 
R   = cond(3); 
Oc  = cond(4); 
Kla = cond(5);
Oh  = cond(6);
Gf  = feed(1);       
Ff  = feed(2);        
Sf  = feed(3);
Ogf = feed(4);


%% Set extracellular metabolite conditions
X  = y(1);
G  = y(2);
F  = y(3);
S  = y(4);
E  = y(5);
Og = y(6);
Ol = y(7);

%% Update bounds and solve for fluxes
[flux,penalty] = solveModel(t,y,INFO);

%% Specific growth rate & uptake rate
mu   = zeros(1,N);
vglc = zeros(1,N);
vfru = zeros(1,N);
vscr = zeros(1,N);
veth = zeros(1,N);
voxy = zeros(1,N);

for i = 1:N
    mu(i)   = flux(i,1);
    vglc(i) = flux(i,2);
    vfru(i) = flux(i,3);
    vscr(i) = flux(i,4);
    veth(i) = flux(i,5);
    voxy(i) = flux(i,6);    
end

%% Calculate model equations
fluxy = [];
dy = []; 

% Saturation gas concentrations  
% (mmol·L-1) <= (mmol·L-1)*(L·Pa·K-1·mmol-1)*K*(atm·Pa-1)*(mol·L-1·atm-1)*(mmol/mol) 
Ols = Og*R*Tr/101325*Oh*1000;

% Extracellular balances·
for i = 1:N
    dX  = -D*X(i)  + mu(i)*X(i);
    dG  = D*(Gf-G) + vglc(i)*X(i);
    dF  = D*(Ff-F) + vfru(i)*X(i);
    dS  = D*(Sf-S) + vscr(i)*X(i);
    dE  = -D*E     + veth(i)*X(i);
    dOg = Dg*(Ogf-Og) - Kla*(Ols-Ol);
    dOl = -D*Ol    + voxy(i)*X(i) + Kla*(Ols-Ol);     
       
    fluxy = [fluxy, mu(i), vglc(i), vfru(i), vscr(i), veth(i), voxy(i)];
    
    dy = [dy dX dG dF dS dE dOg dOl];
    
    dylast = zeros(1, nmodel);
    for j = 1:nmodel
        dylast(i) = penalty(i);
    end
    
    dy = [dy, dylast];
    dy = dy';
end