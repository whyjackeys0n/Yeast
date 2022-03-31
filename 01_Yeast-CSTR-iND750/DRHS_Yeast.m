function [dy,fluxy]= DRHS_Yeast(t, y, INFO)

% Assign values
nmodel = INFO.nmodel;
N = INFO.N;
ns = INFO.ns;
D = INFO.D;
Gf = INFO.Gf;
Ff = INFO.Ff;
Sf = INFO.Sf;

% Set extracellular metabolite conditions
X = y(1);
glc = y(2);
fru = y(3);
scr = y(4);
eth = y(5);

%% Update bounds and solve for fluxes
[flux,penalty] = solveModel(t,y,INFO);
fluxy = [];

%%  
for i = 1:N
    mu(i) = flux(i,1);
    vglc(i) = flux(i,2);
    vfru(i) = flux(i,3);
    vscr(i) = flux(i,4);
    veth(i) = flux(i,5);
    
%     if y(end-nmodel+i) > 1e-3
%         ca_mu(i) = 0;
%         ca_c2(i) = 0;
%         ca_a(i) = 0;
%         ca_e(i) = 0; 
%         ca_bdo(i) = 0;
%         ca_lac(i) = 0;
%         ca_but(i) = 0;
%         penalty(i) = 0;
%     end
    
end

%% Calculate model equations
dy = [];
for i=1:N
    dX = -D*X(i)+mu(i)*X(i);
    dG = D*(Gf-glc)+vglc(i)*X(i);
    dF = D*(Ff-fru)+vfru(i)*X(i);
    dS = D*(Sf-scr)+vscr(i)*X(i);
    dE = -D*eth+veth(i)*X(i);
    
%     dG = D*10+vglc(i)*X(i);
%     dF = D*10+vfru(i)*X(i);
%     dS = D*10+vscr(i)*X(i);
%     dE = -D*10+veth(i)*X(i);
    fluxy = [fluxy mu(i) vglc(i) vfru(i) vscr(i) veth(i)];
    dy = [dy dX dG dF dS dE];

    dylast=zeros(1,nmodel);
    for i=1:nmodel
        dylast(i) = penalty(i);
    end
    dy = [dy dylast];
    dy = dy';
end