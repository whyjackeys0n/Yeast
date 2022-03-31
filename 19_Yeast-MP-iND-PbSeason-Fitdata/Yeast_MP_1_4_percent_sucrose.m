%% Simulation for co-culture seasons data
clear;
clc;

%% Set season number, transferring ratio and initial conditions
seasonnum = 24;
trans = 150;
Xpbi = 250*1e6*4.6e-11*0.8;     % g/L
Xpri = 250*1e6*4.6e-11*0.1;     % g/L
Xchi = 250*1e6*4.6e-11*0.1;     % g/L
Colpb = [255,153,153]./255;
Colpr = [157,195,230]./255;
Colch = [255,217,102]./255;

Yeast_MP_sucrose_function(1,Xpbi,Xpri,Xchi,1);
for i = 1:seasonnum-1
    eval(['load data',num2str(i)]);
    eval(['Yeast_MP_sucrose_function(1,Y',num2str(i),'(end,1)/trans,Y',...
        num2str(i),'(end,2)/trans,Y',num2str(i),'(end,3)/trans,',num2str(i+1),');']);
end
eval(['load data',num2str(i+1)]);

%% Plot 1
for i = 1:seasonnum
%     eval(['figure(',num2str(i),');']);
    eval(['f1=Y',num2str(i),'(:,1)./(Y',num2str(i),'(:,1)+Y',num2str(i),'(:,2)+Y',num2str(i),'(:,3));']);
    eval(['f2=Y',num2str(i),'(:,2)./(Y',num2str(i),'(:,1)+Y',num2str(i),'(:,2)+Y',num2str(i),'(:,3));']);
    eval(['f3=Y',num2str(i),'(:,3)./(Y',num2str(i),'(:,1)+Y',num2str(i),'(:,2)+Y',num2str(i),'(:,3));']);
    eval(['y',num2str(i),'=[f3,f1,f2];']);
%     eval(['area(T',num2str(i),',y',num2str(i),');']);
end

yend = y1(1,:);
figure(1)
for i = 1:seasonnum
    eval(['yend=[yend;y',num2str(i),'(end,:)];']);
end
a1 = area(0:seasonnum,yend);
a1(1).FaceColor = Colch;
a1(2).FaceColor = Colpb;
a1(3).FaceColor = Colpr;
xlim([0,seasonnum]);
ylim([0,1]);
xlabel('Season','FontSize',30);
ylabel('Frequency','FontSize',30);
strains = {'Cheat','Public','Private'};
legend(strains,'FontSize',30)
set(gca,'FontSize',30)
set(gcf,'position',[300,300,1200,1000]);

%% Plot 2
y = [];
T = [];
figure(2)
for i = 1:seasonnum
    eval(['T',num2str(i),'=T',num2str(i),'+24*(',num2str(i),'-1);']);
    eval(['T=[T;T',num2str(i),'];']);
    eval(['y=[y;y',num2str(i),'];']);
end
a2 = area(T,y);
a2(1).FaceColor = Colch;
a2(2).FaceColor = Colpb;
a2(3).FaceColor = Colpr;
xlim([0,T(end)]);
ylim([0,1]);
xlabel('Time [h]','FontSize',30);
ylabel('Frequency','FontSize',30);
strains = {'Cheat','Public','Private'};
legend(strains,'FontSize',30)
set(gca,'FontSize',30)
set(gcf,'position',[300,300,1200,1000]);

%% Plot 3
figure(3)
s1 = Y1(1,1);
s2 = Y1(1,2);
s3 = Y1(1,3);
for i = 1:seasonnum
    eval(['s1=[s1,Y',num2str(i),'(end,1)];']);
end
for i = 1:seasonnum
    eval(['s2=[s2,Y',num2str(i),'(end,2)];']);
end
for i = 1:seasonnum
    eval(['s3=[s3,Y',num2str(i),'(end,3)];']);
end
s = s1 + s2+ s3;
semilogy(0:seasonnum,s,'-ko','LineWidth',5,'MarkerFaceColor','k','MarkerSize',15);
xlim([0,seasonnum]);
xlabel('Season','FontSize',30);
ylabel('Concentration [g/L]','FontSize',30);
set(gca,'FontSize',30);
set(gcf,'position',[300,300,1200,1000]);

%% Plot 4
figure(4)
scr = [];
for i = 1:seasonnum
    eval(['sc=Y',num2str(i),'(:,10);']);
    scr=[scr;sc];
end
plot(T,scr,'-k','LineWidth',3);
xlabel('Time [h]','FontSize',30);
ylabel('Invertase [g/L]','FontSize',30);
set(gca,'FontSize',30);
set(gcf,'position',[300,300,1200,1000]);

%% Plot 5
figure(5)
biomasspb = [];
biomasspr = [];
biomassch = [];
for i = 1:seasonnum
    eval(['biopb=Y',num2str(i),'(:,1);']);
    eval(['biopr=Y',num2str(i),'(:,2);']);
    eval(['bioch=Y',num2str(i),'(:,3);']);
    biomasspb=[biomasspb;biopb];
    biomasspr=[biomasspr;biopr];
    biomassch=[biomassch;bioch];
end
plot(T,biomasspb,'-','LineWidth',3);
hold on
plot(T,biomasspr,'-','LineWidth',3);
hold on
plot(T,biomassch,'-','LineWidth',3);
species = {'Public','Private','Cheat'};
legend(species,'FontSize',30)
xlabel('Time [h]','FontSize',30);
ylabel('Biomass [g/L]','FontSize',30);
set(gca,'FontSize',30);
set(gcf,'position',[300,300,1200,1000]);

%% Plot settings
% species = {'0.0625% Sucrose','0.25% Sucrose'};
% legend(species,'FontSize',30)
% fh = figure(1);
% fh.WindowState = 'maximized';
