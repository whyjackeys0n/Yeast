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

function [dy, fluxy, coeff, ca_mu, er_mu, ca_c, ca_c2, ca_a, ca_e, ca_bdo, ca_lac, ca_h2, er_g, er_c2, er_but, er_a, er_lac, er_h2]= DRHS(t, y, INFO)
% Assign values
nmodel = INFO.nmodel;
N = INFO.N;
ns = INFO.ns;
condit = INFO.condit;

% Set conditions
klac = condit(1);
Hc = condit(2);
pc = condit(3);
pc2 = condit(4);
klac2 = condit(5);
Hc2 = condit(6);
tr = condit(7);
zs = condit(8);
% N = condit(9);
% ns = condit(10);
ug0 = condit(11);
ul = condit(12);
dl = condit(13);
eg = condit(14);
el = condit(15);
% cgi = condit(16);
pL = condit(17);
% c2gi = condit(18);
Area = condit(19);
Qmedia = condit(20);
D = condit(21);
densityL = condit(22);
g = condit(23);
viscosityL = condit(24);
db0 = condit(25);
pn = condit(26);
klc = condit(27);
surfaceTL = condit(29);
gi = condit(30);
ph2 = condit(31);
Hh2 = condit(32);
klc2 = klc;

% Define extracellular state variables
for i=1:N
    cg(i) = y(1+(i-1)*ns);  
    c2g(i) = y(2+(i-1)*ns);  
    h2g(i) = y(3+(i-1)*ns); 
    cl(i) = y(4+(i-1)*ns);
    c2l(i) = y(5+(i-1)*ns);
    h2l(i) = y(6+(i-1)*ns);
    p(i) = y(7+(i-1)*ns);  
    jg(i) = y(8+(i-1)*ns);  
    db(i) = y(9+(i-1)*ns); 
    ub(i) = y(10+(i-1)*ns); 
    eg(i) = y(11+(i-1)*ns); 
end

Xca = y(N*ns+1);
Xer = y(N*ns+2);
A = y(N*ns+3);
E = y(N*ns+4);
BDO = y(N*ns+5);
Lac = y(N*ns+6);
G = y(N*ns+7);
BUT = y(N*ns+8);

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
    ca_h2(i) = flux(i,8);
    
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
    er_h2(i) = flux(N+i,7);  
    
    if y(end-nmodel+N+i) > 1e-3
        er_mu(i) = 0;
        er_c2(i) = 0;
        er_but(i) = 0;
        er_lac(i) = 0;
        er_h2(i) = 0;
        penalty(N+i) = 0;        
    end
end
aca_mu = sum(ca_mu)/N;
aca_a = sum(ca_a)/N;
aca_e = sum(ca_e)/N;
aca_bdo = sum(ca_bdo)/N;
aca_lac = sum(ca_lac)/N;
aer_mu = sum(er_mu)/N;
aer_g = sum(er_g)/N;
aer_but = sum(er_but)/N;
aer_a = sum(er_a)/N;
aer_lac = sum(er_lac)/N;
%% Dynamics

dy = [];

i=1:N;
% Saturation gas concentrations
    cls(i) = cg(i)*8.314*tr*Hc/1.013e5*1000.*p(i)/1.013e5;
    c2ls(i) = c2g(i)*8.314*tr*Hc2/1.013e5*1000.*p(i)/1.013e5;
    h2ls(i) = h2g(i)*8.314*tr*Hh2/1.013e5*1000.*p(i)/1.013e5;
    n2g(i) = p(i)*0.3/8.314/tr;
    el(i) = 1-eg(i);
    r(i) = (cg(i)+c2g(i)+h2g(i)+n2g(i)).*jg(i)/(cg(1)+c2g(1)+h2g(1)+n2g(1))/jg(1);
        
i=1;
        p_bot = p(1);
        cgi = pc*p_bot/8.314/tr;  
        c2gi = pc2*p_bot/8.314/tr; 
        h2gi = ph2*p_bot/8.314/tr;
        n2gi = pn*p_bot/8.314/tr;
        dbi = db0*nthroot((pL/p_bot),3);
        jgi = ug0*pL/p_bot;
        ubi = 0.33*(g^0.76)*((densityL/viscosityL)^0.52)*(((dbi/2)*1e-3).^1.28);
        egi = jgi/ubi;
        ngi = (cgi+c2gi+h2gi+n2gi)*jgi;
        
        dp(i) = -densityL*g*el(i) - (p(i+1) - p(i))/zs;
        ddb(i) = (pL/p(i))*((cg(i)+c2g(i)+n2g(i))*jg(i)/ngi) - (db(i)/db0)^3;
        djg(i) = ug0*(pL/p(i))*r(i) - jg(i);
        dub(i) = 0.33*(g^0.76)*((densityL/viscosityL)^0.52)*(((db(i)/2)*1e-3).^1.28) - ub(i);
        deg(i) = (sqrt((1.05*((jg(i)*3600)+ul)+(ub(i)*3600))^2-4*(ub(i)*3600)*(jg(i)*3600)) + (-1.05*((jg(i)*3600)+ul)-(ub(i)*3600)))/(-2*(ub(i)*3600)) - eg(i);
        klac(i) = 6*eg(i)/(1-eg(i))/(db(i)/1000)*klc;
        klac2(i) = 6*eg(i)/(1-eg(i))/(db(i)/1000)*klc2;
        
        jgd(i) = ((jg(i) - jgi)*3600)/zs;
        eli = 1 - egi;
        eld(i) = (el(i) - eli)/zs;
        egd(i) = (eg(i) - egi)/zs;
        cgd(i) = (cg(i)-cgi)/zs;
        c2gd(i) = (c2g(i)-c2gi)/zs;  
        h2gd(i) = (h2g(i)-h2gi)/zs; 
        cld(i) = 0;
        c2ld(i) = 0;
        h2ld(i) = 0;
        cld2(i) = (cl(i+1)-cl(i))/zs^2;
        c2ld2(i) = (c2l(i+1)-c2l(i))/zs^2;
        h2ld2(i) = (h2l(i+1)-h2l(i))/zs^2;

        dcg(i) = -(jg(i)*3600)*cgd(i)/eg(i)-cg(i)*jgd(i)/eg(i)-klac(i)/eg(i)*(cls(i)-cl(i));
        dc2g(i) = -(jg(i)*3600)*c2gd(i)/eg(i)-c2g(i)*jgd(i)/eg(i)-klac2(i)/eg(i)*(c2ls(i)-c2l(i));
        dh2g(i) = -(jg(i)*3600)*h2gd(i)/eg(i)-h2g(i)*jgd(i)/eg(i)-klac(i)/eg(i)*(h2ls(i)-h2l(i));
        dcl(i) = ca_c(i)*Xca+klac(i)*(cls(i)-cl(i))/el(i)+dl*cld2(i)+dl*eld(i)*cld(i)/el(i)-(ul-Qmedia/Area)*cld(i)/el(i)-(ul-Qmedia/Area)*cl(i)*eld(i)/el(i); 
        dc2l(i) = ca_c2(i)*Xca+er_c2(i)*Xer+klac2(i)*(c2ls(i)-c2l(i))/el(i)+dl*c2ld2(i)+dl*c2ld(i)*eld(i)/el(i)-(ul-Qmedia/Area)*c2ld(i)/el(i)-(ul-Qmedia/Area)*c2l(i)*eld(i)/el(i); 
        dh2l(i) = ca_h2(i)*Xca+er_h2(i)*Xer+klac(i)*(h2ls(i)-h2l(i))/el(i)+dl*h2ld2(i)+dl*h2ld(i)*eld(i)/el(i)-(ul-Qmedia/Area)*h2ld(i)/el(i)-(ul-Qmedia/Area)*h2l(i)*eld(i)/el(i); 

i=2:N-1;
        dp(i) = -densityL*g*el(i) - (p(i+1) - p(i))/zs;
        ddb(i) = (pL./p(i)).*((cg(i)+c2g(i)+n2g(i)).*jg(i)./ngi) - (db(i)./db0).^3;
        djg(i) = ug0*(pL./p(i)).*r(i) - jg(i);
        dub(i) = 0.33*(g^0.76)*((densityL/viscosityL)^0.52)*(((db(i)/2)*1e-3).^1.28) - ub(i);
        deg(i) = (sqrt((1.05*((jg(i)*3600)+ul)+(ub(i)*3600)).^2-4*(ub(i)*3600).*(jg(i)*3600)) + (-1.05*((jg(i)*3600)+ul)-(ub(i)*3600)))./(-2*(ub(i)*3600)) - eg(i);
        klac(i) = 6*eg(i)./(1-eg(i))./(db(i)/1000)*klc;
        klac2(i) = 6*eg(i)./(1-eg(i))./(db(i)/1000)*klc2;
        
        jgd(i) = ((jg(i) - jg(i-1))*3600)/zs;
        eld(i) = (el(i) - el(i-1))/zs;
        egd(i) = (eg(i) - eg(i-1))/zs;
        cgd(i) = (cg(i)-cg(i-1))/zs;
        c2gd(i) = (c2g(i)-c2g(i-1))/zs;
        h2gd(i) = (h2g(i)-h2g(i-1))/zs;
        cld(i) = (cl(i+1)-cl(i))/zs;
        c2ld(i) = (c2l(i+1)-c2l(i))/zs;
        h2ld(i) = (h2l(i+1)-h2l(i))/zs;
        cld2(i) = (cl(i+1)-2*cl(i)+cl(i-1))/zs^2;
        c2ld2(i) = (c2l(i+1)-2*c2l(i)+c2l(i-1))/zs^2;
        h2ld2(i) = (h2l(i+1)-2*h2l(i)+h2l(i-1))/zs^2;
        
        dcg(i) = -(jg(i)*3600).*cgd(i)./eg(i)-cg(i).*jgd(i)./eg(i)-klac(i)./eg(i).*(cls(i)-cl(i));
        dc2g(i) = -(jg(i)*3600).*c2gd(i)./eg(i)-c2g(i).*jgd(i)./eg(i)-klac2(i)./eg(i).*(c2ls(i)-c2l(i));
        dh2g(i) = -(jg(i)*3600).*h2gd(i)./eg(i)-h2g(i).*jgd(i)./eg(i)-klac(i)./eg(i).*(h2ls(i)-h2l(i)); 
        dcl(i) = ca_c(i)*Xca+klac(i).*(cls(i)-cl(i))./el(i)+dl*cld2(i)+dl*eld(i).*cld(i)./el(i)-(ul-Qmedia/Area)*cld(i)./el(i)-(ul-Qmedia/Area)*cl(i).*eld(i)./el(i);  
        dc2l(i) = ca_c2(i)*Xca+er_c2(i)*Xer+klac2(i).*(c2ls(i)-c2l(i))./el(i)+dl*c2ld2(i)+dl*c2ld(i).*eld(i)./el(i)-(ul-Qmedia/Area)*c2ld(i)./el(i)-(ul-Qmedia/Area)*c2l(i).*eld(i)./el(i); 
        dh2l(i) = ca_h2(i)*Xca+er_h2(i)*Xer+klac(i).*(h2ls(i)-h2l(i))./el(i)+dl*h2ld2(i)+dl*h2ld(i).*eld(i)./el(i)-(ul-Qmedia/Area)*h2ld(i)./el(i)-(ul-Qmedia/Area)*h2l(i).*eld(i)./el(i); 

        
i=N;
        dp(i) = pL - p(i);
        ddb(i) = (pL/p(i))*((cg(i)+c2g(i)+n2g(i))*jg(i)/ngi) - (db(i)/db0)^3;
        djg(i) = ug0*(pL/p(i))*r(i) - jg(i);
        dub(i) = 0.33*(g^0.76)*((densityL/viscosityL)^0.52)*(((db(i)/2)*1e-3).^1.28) - ub(i);
        deg(i) = (sqrt((1.05*((jg(i)*3600)+ul)+(ub(i)*3600))^2-4*(ub(i)*3600)*(jg(i)*3600)) + (-1.05*((jg(i)*3600)+ul)-(ub(i)*3600)))/(-2*(ub(i)*3600)) - eg(i);
        klac(i) = 6*eg(i)/(1-eg(i))/(db(i)/1000)*klc;
        klac2(i) = 6*eg(i)/(1-eg(i))/(db(i)/1000)*klc2;

        jgd(i) = ((jg(i) - jg(i-1))*3600)/zs;
        eld(i) = (el(i) - el(i-1))/zs;
        egd(i) = (eg(i) - eg(i-1))/zs;
        cgd(i) = (cg(i)-cg(i-1))/zs;            
        c2gd(i) = (c2g(i)-c2g(i-1))/zs;
        h2gd(i) = (h2g(i)-h2g(i-1))/zs;
        c2lr = c2l(1)*(abs(ul)*Area-Qmedia)/(abs(ul)*Area); 
        clr = 0;
        h2lr = 0;
        cln = ((ul-Qmedia/Area)*zs*clr/dl-cl(i))/((ul-Qmedia/Area)*zs/dl-1); 
        c2ln = ((ul-Qmedia/Area)*zs*c2lr/dl-c2l(i))/((ul-Qmedia/Area)*zs/dl-1);
        h2ln = ((ul-Qmedia/Area)*zs*h2lr/dl-h2l(i))/((ul-Qmedia/Area)*zs/dl-1);
        coeff=[cln, c2ln, h2ln];
        cld(i) = (cln-cl(i))/zs;
        c2ld(i) = (c2ln-c2l(i))/zs;
        h2ld(i) = (h2ln-h2l(i))/zs;
        cld2(i) = (cln-2*cl(i)+cl(i-1))/zs^2;
        c2ld2(i) = (c2ln-2*c2l(i)+c2l(i-1))/zs^2;
        h2ld2(i) = (h2ln-2*h2l(i)+h2l(i-1))/zs^2;

        dcg(i) = -(jg(i)*3600)*cgd(i)/eg(i)-cg(i)*jgd(i)/eg(i)-klac(i)/eg(i)*(cls(i)-cl(i));
        dc2g(i) = -(jg(i)*3600)*c2gd(i)/eg(i)-c2g(i)*jgd(i)/eg(i)-klac2(i)/eg(i)*(c2ls(i)-c2l(i));
        dh2g(i) = -(jg(i)*3600)*h2gd(i)/eg(i)-h2g(i)*jgd(i)/eg(i)-klac(i)/eg(i)*(h2ls(i)-h2l(i));
        dcl(i) = ca_c(i)*Xca+klac(i)*(cls(i)-cl(i))/el(i)+dl*cld2(i)+dl*eld(i)*cld(i)/el(i)-(ul-Qmedia/Area)*cld(i)/el(i)-(ul-Qmedia/Area)*cl(i)*eld(i)/el(i);  
        dc2l(i) = ca_c2(i)*Xca+er_c2(i)*Xer+klac2(i)*(c2ls(i)-c2l(i))/el(i)+dl*c2ld2(i)+dl*c2ld(i)*eld(i)/el(i)-(ul-Qmedia/Area)*c2ld(i)/el(i)-(ul-Qmedia/Area)*c2l(i)*eld(i)/el(i);   
        dh2l(i) = ca_h2(i)*Xca+er_h2(i)*Xer+klac(i)*(h2ls(i)-h2l(i))/el(i)+dl*h2ld2(i)+dl*h2ld(i)*eld(i)/el(i)-(ul-Qmedia/Area)*h2ld(i)/el(i)-(ul-Qmedia/Area)*h2l(i)*eld(i)/el(i); 

for i=1:N
    dyTmp = [dcg(i); dc2g(i); dh2g(i); dcl(i); dc2l(i); dh2l(i); dp(i); djg(i); ddb(i); dub(i); deg(i)];
    dyTmp1 = reshape(dyTmp,[],1);
    dy=[dy dyTmp1'];
end

dXca = aca_mu*Xca - D*Xca;
dXer = aer_mu*Xer - D*Xer;           
dA =  aca_a*Xca + aer_a*Xer - D*A;               
dE =  aca_e*Xca - D*E; 
dBDO =  aca_bdo*Xca - D*BDO; 
dL =  aca_lac*Xca + aer_lac*Xer - D*Lac;
dG = aer_g*Xer + D*gi - D*G;
dBUT = aer_but*Xer - D*BUT;      

dy = [dy dXca dXer dA dE dBDO dL dG dBUT];

dylast=zeros(1,nmodel);
for i=1:nmodel
    dylast(i) = penalty(i);
end
dy = [dy dylast];
dy = dy';
end