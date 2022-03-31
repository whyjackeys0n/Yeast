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
Xpr = y(1);
Xch = y(2);
G   = y(3);
F   = y(4);
S   = y(5);
E   = y(6);
Og  = y(7);
Ol  = y(8);

%% Update bounds and solve for fluxes
[flux,penalty] = solveModel(t,y,INFO);

%% Specific growth rate & uptake rate

for i = 1:N
    pr_mu(i)   = flux(i,1);
    pr_vglc(i) = flux(i,2);
    pr_vfru(i) = flux(i,3);
    pr_vscr(i) = flux(i,4);
    pr_veth(i) = flux(i,5);
    pr_voxy(i) = flux(i,6);    
    
    ch_mu(i)   = flux(N+i,1);
    ch_vglc(i) = flux(N+i,2);
    ch_vfru(i) = flux(N+i,3);
    ch_vscr(i) = flux(N+i,4);
    ch_veth(i) = flux(N+i,5);
    ch_voxy(i) = flux(N+i,6);
end

%% Calculate model equations
fluxy = [];
dy = []; 

% Saturation gas concentrations  
% (mmol·L-1) <= (mmol·L-1)*(L·Pa·K-1·mmol-1)*K*(atm·Pa-1)*(mol·L-1·atm-1)*(mmol/mol) 
Ols = Og*R*Tr/101325*Oh*1000;

% Extracellular balances·
for i = 1:N
    dXpr = -D*Xpr(i) + pr_mu(i)*Xpr(i);
    dXch = -D*Xch(i) + ch_mu(i)*Xch(i);
    dG   = D*(Gf-G)  + pr_vglc(i)*Xpr(i) + ch_vglc(i)*Xch(i);
    dF   = D*(Ff-F)  + pr_vfru(i)*Xpr(i) + ch_vfru(i)*Xch(i);
    dS   = D*(Sf-S)  + pr_vscr(i)*Xpr(i) + ch_vscr(i)*Xch(i);
    dE   = -D*E      + pr_veth(i)*Xpr(i) + ch_veth(i)*Xch(i);
    dOg  = Dg*(Ogf-Og) - Kla*(Ols-Ol);
    dOl  = -D*Ol     + pr_voxy(i)*Xpr(i) + ch_voxy(i)*Xch(i) + Kla*(Ols-Ol);    
       
    fluxy = [fluxy, pr_mu(i), pr_vglc(i), pr_vfru(i), pr_vscr(i), pr_veth(i), pr_voxy(i)...
                  , ch_mu(i), ch_vglc(i), ch_vfru(i), ch_vscr(i), ch_veth(i), ch_voxy(i)];
    
    dy = [dy dXpr dXch dG dF dS dE dOg dOl];
    
    dylast = zeros(1, nmodel);
    for j = 1:nmodel
        dylast(i) = penalty(i);
    end
    
    dy = [dy, dylast];
    dy = dy';
end