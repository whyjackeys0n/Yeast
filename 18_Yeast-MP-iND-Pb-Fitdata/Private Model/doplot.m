load data1.mat
load data2.mat

FS = 30;

figure(1)
subplot(2,4,1)
plot(T1,Y1(:,2),'linewidth',2.5);
hold on
plot(T2,Y2(:,2),'linewidth',2.5);
ylabel('Biomass (g/L)');
set(gca,'FontSize',FS)

subplot(2,4,2)
plot(T1,Y1(:,6),'linewidth',2.5);
hold on
plot(T2,Y2(:,6),'linewidth',2.5);
ylabel('Sucrose (mmol/L)');
set(gca,'FontSize',FS)

subplot(2,4,3)
plot(T1,Y1(:,4),'linewidth',2.5);
hold on
plot(T2,Y2(:,4),'linewidth',2.5);
ylabel('Glucose (mmol/L)');
set(gca,'FontSize',FS)

subplot(2,4,4)
plot(T1,Y1(:,5),'linewidth',2.5);
hold on
plot(T2,Y2(:,5),'linewidth',2.5);
ylabel('Fructose (mmol/L)');
set(gca,'FontSize',FS)

subplot(2,4,5)
plot(T1,Y1(:,7),'linewidth',2.5);
hold on
plot(T2,Y2(:,7),'linewidth',2.5);
ylabel('Ethanol (mmol/L)');
set(gca,'FontSize',FS)

subplot(2,4,6)
plot(T1,Y1(:,9),'linewidth',2.5);
hold on
plot(T2,Y2(:,9),'linewidth',2.5);
ylabel('Dissolved O2 (mmol/L)');
set(gca,'FontSize',FS)

subplot(2,4,7)
plot(T1,Y1(:,10),'linewidth',2.5);
hold on
plot(T2,Y2(:,10),'linewidth',2.5);
ylabel('Invertase (g/L)');
set(gca,'FontSize',FS)

figure(2)
subplot(2,3,1)
plot(T1,flux1(:,7),'linewidth',2.5);
hold on
plot(T2,flux2(:,7),'linewidth',2.5);
ylabel('Growth (h-1)');
set(gca,'FontSize',FS)

subplot(2,3,2)
plot(T1,flux1(:,10),'linewidth',2.5);
hold on
plot(T2,flux2(:,10),'linewidth',2.5);
ylabel('Sucrose (mmol/g/h)');
set(gca,'FontSize',FS)

subplot(2,3,3)
plot(T1,flux1(:,8),'linewidth',2.5);
hold on
plot(T2,flux2(:,8),'linewidth',2.5);
ylabel('Glucose (mmol/g/h)');
set(gca,'FontSize',FS)

subplot(2,3,4)
plot(T1,flux1(:,9),'linewidth',2.5);
hold on
plot(T2,flux2(:,9),'linewidth',2.5);
ylabel('Fructose (mmol/g/h)');
set(gca,'FontSize',FS)

subplot(2,3,5)
plot(T1,flux1(:,11),'linewidth',2.5);
hold on
plot(T2,flux2(:,11),'linewidth',2.5);
ylabel('Ethanol (mmol/g/h)');
set(gca,'FontSize',FS)

subplot(2,3,6)
plot(T1,flux1(:,12),'linewidth',2.5);
hold on
plot(T2,flux2(:,12),'linewidth',2.5);
ylabel('Dissolved O2 (mmol/g/h)');
set(gca,'FontSize',FS)
