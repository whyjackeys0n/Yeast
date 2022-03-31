function [dy,fluxy]= DRHS(t, y, INFO)

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
gi = condit(13);
Di = condit(14);
Dgi = condit(15); 
D = Di;
Dg = Dgi;

Xca = y(1);
Xer = y(2);
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
G = y(13);
BUT = y(14);

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
    
    if y(end-nmodel+i) > 1e-3
        ca_mu(i) = 0;
        ca_c2(i) = 0;
        ca_a(i) = 0;
        ca_e(i) = 0; 
        ca_bdo(i) = 0;
        ca_lac(i) = 0;
        penalty(i) = 0;
    end

    er_mu(i) = flux(N+i,1);
    er_g(i) = flux(N+i,2);   
    er_c2(i) = flux(N+i,3);   
    er_but(i) = flux(N+i,4);   
    er_a(i) = flux(N+i,5);     
    er_lac(i) = flux(N+i,6);
    er_h(i) = flux(N+i,7);  
    
    if y(end-nmodel+N+i) > 1e-3
        er_mu(i) = 0;
        er_c2(i) = 0;
        er_but(i) = 0;
        er_lac(i) = 0;
        er_h(i) = 0;
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
dXer = er_mu(i)*Xer - D*Xer;
dcg = Dg*(cgi-cg)-klac*(cls(i)-cl(i));
dc2g = Dg*(c2gi-c2g)-klac2*(c2ls(i)-c2l(i));
dhg = Dg*(hgi-hg)-klah*(hls(i)-hl(i));
dcl(i) = ca_c(i)*Xca + klac*(cls(i)-cl(i)) - D*cl(i);
dc2l(i) = ca_c2(i)*Xca +er_c2(i)*Xer + klac2*(c2ls(i)-c2l(i)) - D*c2l(i);
dhl(i) = ca_h(i)*Xca +er_h(i)*Xer + klah*(hls(i)-hl(i)) - D*hl(i);
dA =  ca_a(i)*Xca + er_a(i)*Xer - D*A;               
dE =  ca_e(i)*Xca - D*E; 
dBDO =  ca_bdo(i)*Xca - D*BDO; 
dL =  ca_lac(i)*Xca + er_lac(i)*Xer - D*Lac;
dG = er_g(i)*Xer + D*gi - D*G;
dBUT = er_but(i)*Xer - D*BUT;

fluxy = [fluxy ca_mu(i) ca_c(i) ca_c2(i) ca_a(i) ca_e(i) ca_bdo(i) ca_lac(i) ca_h(i) er_mu(i) er_g(i) er_c2(i) er_but(i) er_a(i) er_lac(i) er_h(i)];
dy = [dXca dXer dcg dc2g dhg dcl dc2l dhl dA dE dBDO dL dG dBUT];

dylast=zeros(1,nmodel);
for i=1:nmodel
    dylast(i) = penalty(i);
end
dy = [dy dylast];
dy = dy';
end