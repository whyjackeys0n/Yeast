function [dy,fluxy]= DRHS_Ckluyveri(t, y, INFO)

% Assign values
nmodel = INFO.nmodel;
N = INFO.N;
ns = INFO.ns;
condit = INFO.condit;

% Set conditions
klac = condit(1);
Hc = condit(2);
klac2 = condit(3);
Hc2 = condit(4);
klah = condit(5);
Hh = condit(6);
tr = condit(7);
eg = condit(8);
el = condit(9);
cgi = condit(10);
c2gi = condit(11);
hgi = condit(12);
D = condit(13);
suci = condit(14);
proi = condit(15);
croi = condit(16);
vinylai = condit(17);
propii = condit(18);
buti = condit(19);
Dg = condit(20);

% Define extracellular state variables
Xca = y(1);
Xck = y(2);
cg = y(3); 
c2g = y(4);
hg = y(5);
cl = y(6);
c2l = y(7);
hl = y(8);
A = y(9);
E = y(10);
BDO = y(11);
Lac = y(12);
suc = y(13);
pro = y(14);
cro = y(15);
vinyla = y(16);
propi = y(17);
but = y(18);
pen = y(19);
hex = y(20);

%% Update bounds and solve for fluxes
[flux,penalty] = solveModel(t,y,INFO);
fluxy = [];
%%
   
for i = 1:N

    ca_mu(i) = flux(i,1); 
    ca_c(i) = flux(i,2);     
    ca_c2(i) = flux(i,3);   
    ca_a(i) = flux(i,4);      
    ca_e(i) = flux(i,5);      
    ca_bdo(i) = flux(i,6);
    ca_lac(i) = flux(i,7);
    ca_h(i) = flux(i,8);
    ca_but(i) = flux(i,9);
    
    if y(end-nmodel+i) > 1e-3
        ca_mu(i) = 0;
        ca_c2(i) = 0;
        ca_a(i) = 0;
        ca_e(i) = 0; 
        ca_bdo(i) = 0;
        ca_lac(i) = 0;
        ca_but(i) = 0;
        penalty(i) = 0;
    end

    ck_mu(i) = flux(N+i,1);
    ck_a(i) = flux(N+i,2);   
    ck_e(i) = flux(N+i,3);   
    ck_suc(i) = flux(N+i,4);   
    ck_pro(i) = flux(N+i,5);     
    ck_cro(i) = flux(N+i,6);
    ck_vinyla(i) = flux(N+i,7);  
    ck_propi(i) = flux(N+i,8); 
    ck_but(i) = flux(N+i,9); 
    ck_pen(i) = flux(N+i,10);
    ck_hex(i) = flux(N+i,11);
    ck_h2(i) = flux(N+i,12);
    
    if y(end-nmodel+N+i) > 1e-3
        ck_mu(i) = 0;
        ck_but(i) = 0;
        penalty(N+i) = 0;        
    end
end

%% Dynamics

dy = [];

i=1:N;
% Saturation gas concentrations
cls(i) = cg*8.314*tr*Hc/1.013e5*1000;
c2ls(i) = c2g*8.314*tr*Hc2/1.013e5*1000;
hls(i) = hg*8.314*tr*Hc2/1.013e5*1000;

% Extracellular balances     
i=1:N; 
dXca = ca_mu(i)*Xca - D*Xca;
dXck = ck_mu(i)*Xck - D*Xck;
dcg = Dg*(cgi-cg)-klac*(cls(i)-cl(i));
dc2g = Dg*(c2gi-c2g)-klac2*(c2ls(i)-c2l(i));
dhg = Dg*(hgi-hg)-klah*(hls(i)-hl(i));
dcl(i) = ca_c(i)*Xca + klac*(cls(i)-cl(i)) - D*cl(i);
dc2l(i) = ca_c2(i)*Xca + klac2*(c2ls(i)-c2l(i)) - D*c2l(i);
dhl(i) = ck_h2(i)*Xck + ca_h(i)*Xca + klah*(hls(i)-hl(i)) - D*hl(i);
dA =  ca_a(i)*Xca + ck_a(i)*Xck - D*A;               
dE =  ca_e(i)*Xca + ck_e(i)*Xck - D*E; 
dBDO =  ca_bdo(i)*Xca - D*BDO; 
dL =  ca_lac(i)*Xca - D*Lac;
dsuc = ck_suc(i)*Xck + D*suci - D*suc;
dpro = ck_pro(i)*Xck + D*proi - D*pro;
dcro = ck_cro(i)*Xck + D*croi - D*cro;
dvinyla = ck_vinyla(i)*Xck + D*vinylai - D*vinyla;
dpropi = ck_propi(i)*Xck + D*propii - D*propi;
dbut = ck_but(i)*Xck + ca_but(i)*Xca + D*buti - D*but;
dpen = ck_pen(i)*Xck - D*pen;
dhex = ck_hex(i)*Xck - D*hex;

fluxy = [fluxy ca_mu(i) ca_c(i) ca_c2(i) ca_a(i) ca_e(i) ca_bdo(i) ca_lac(i) ca_h(i) ca_but(i) ck_mu(i) ck_a(i) ck_e(i) ck_suc(i) ck_pro(i) ck_cro(i) ck_vinyla(i) ck_propi(i) ck_but(i) ck_pen(i) ck_hex(i) ck_h2(i)];
dy = [dXca dXck dcg dc2g dhg dcl dc2l dhl dA dE dBDO dL dsuc dpro dcro dvinyla dpropi dbut dpen dhex];

dylast=zeros(1,nmodel);
for i=1:nmodel
    dylast(i) = penalty(i);
end
dy = [dy dylast];
dy = dy';
end