function [dy,fluxy]= DRHS_Yeast(t, y, INFO)

%% Assign values
nmodel = INFO.nmodel;
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
Xpb = y(1);
Xpr = y(2);
Xch = y(3);
G   = y(4);
F   = y(5);
S   = y(6);
E   = y(7);
Og  = y(8);
Ol  = y(9);
Inv = y(10);

%% Update bounds and solve for fluxes
[flux,penalty] = solveModel(t,y,INFO);

%% Specific growth rate & uptake rate
pb_mu   = flux(1,1);
pb_vglc = flux(1,2);
pb_vfru = flux(1,3);
pb_vscr = flux(1,4);
pb_veth = flux(1,5);
pb_voxy = flux(1,6);

pr_mu   = flux(2,1);
pr_vglc = flux(2,2);
pr_vfru = flux(2,3);
pr_vscr = flux(2,4);
pr_veth = flux(2,5);
pr_voxy = flux(2,6);

ch_mu   = flux(3,1);
ch_vglc = flux(3,2);
ch_vfru = flux(3,3);
ch_vscr = flux(3,4);
ch_veth = flux(3,5);
ch_voxy = flux(3,6);

%% Calculate model equations
fluxy = [];
dy = []; 

% Saturation gas concentrations  
% (mmol·L-1) <= (mmol·L-1)*(L·Pa·K-1·mmol-1)*K/(Pa·atm-1)*(mol·L-1·atm-1)*(mmol/mol) 
Ols = Og*R*Tr/101325*Oh*1000;

% Extracellular balances
dXpb = -D*Xpb + pb_mu*Xpb;
dXpr = -D*Xpr + pr_mu*Xpr;
dXch = -D*Xch + ch_mu*Xch;
dG   = D*(Gf-G)  + pb_vglc*Xpb + pr_vglc*Xpr + ch_vglc*Xch + Inv*2;
dF   = D*(Ff-F)  + pb_vfru*Xpb + pr_vfru*Xpr + ch_vfru*Xch;
dS   = D*(Sf-S)  + pb_vscr*Xpb + pr_vscr*Xpr + ch_vscr*Xch - Inv;
dE   = -D*E      + pb_veth*Xpb + pr_veth*Xpr + ch_veth*Xch;
dOg  = Dg*(Ogf-Og) - Kla*(Ols-Ol);
dOl  = -D*Ol     + pb_voxy*Xpb + pr_voxy*Xpr + ch_voxy*Xch + Kla*(Ols-Ol);    
dInv = 0.05*Xpb + 0.005*G/(0.02+G)/(1+G/0.07)*Xpb - 0.3*Inv;

fluxy = [fluxy, pb_mu, pb_vglc, pb_vfru, pb_vscr, pb_veth, pb_voxy...
                pr_mu, pr_vglc, pr_vfru, pr_vscr, pr_veth, pr_voxy...
              , ch_mu, ch_vglc, ch_vfru, ch_vscr, ch_veth, ch_voxy];

dy = [dy dXpb dXpr dXch dG dF dS dE dOg dOl dInv];

dylast = zeros(1, nmodel);
for j = 1:nmodel
    dylast(1) = penalty(1);
end

dy = [dy, dylast];
dy = dy';
