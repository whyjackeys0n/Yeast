load data1.mat
load data2.mat

FS = 20;

figure(1)
subplot(2,3,1)
plot(T1,Y1(:,2),'linewidth',2.5);
hold on
plot(T2,Y2(:,2),'linewidth',2.5);
ylabel('Biomass (g/L)');
set(gca,'FontSize',FS)

subplot(2,3,2)
plot(T1,Y1(:,6),'linewidth',2.5);
hold on
plot(T2,Y2(:,6),'linewidth',2.5);
ylabel('Sucrose (mmol/L)');
set(gca,'FontSize',FS)

subplot(2,3,3)
plot(T1,Y1(:,7),'linewidth',2.5);
hold on
plot(T2,Y2(:,7),'linewidth',2.5);
ylabel('Ethanol (mmol/L)');
set(gca,'FontSize',FS)

subplot(2,3,4)
plot(T1,Y1(:,9),'linewidth',2.5);
hold on
plot(T2,Y2(:,9),'linewidth',2.5);
ylabel('Dissolved O2 (mmol/L)');
set(gca,'FontSize',FS)

subplot(2,3,5)
plot(T1,Y1(:,10),'linewidth',2.5);
hold on
plot(T2,Y2(:,10),'linewidth',2.5);
ylabel('Ammonia (mmol/L)');
set(gca,'FontSize',FS)

figure(2)
subplot(2,3,1)
plot(T1,flux1(:,8),'linewidth',2.5);
hold on
plot(T2,flux2(:,8),'linewidth',2.5);
ylabel('Growth (h-1)');
set(gca,'FontSize',FS)

subplot(2,3,2)
plot(T1,flux1(:,11),'linewidth',2.5);
hold on
plot(T2,flux2(:,11),'linewidth',2.5);
ylabel('Sucrose (mmol/g/h)');
set(gca,'FontSize',FS)

subplot(2,3,3)
plot(T1,flux1(:,12),'linewidth',2.5);
hold on
plot(T2,flux2(:,12),'linewidth',2.5);
ylabel('Ethanol (mmol/g/h)');
set(gca,'FontSize',FS)

subplot(2,3,4)
plot(T1,flux1(:,13),'linewidth',2.5);
hold on
plot(T2,flux2(:,13),'linewidth',2.5);
ylabel('Dissolved O2 (mmol/g/h)');
set(gca,'FontSize',FS)

subplot(2,3,5)
plot(T1,flux1(:,14),'linewidth',2.5);
hold on
plot(T2,flux2(:,14),'linewidth',2.5);
ylabel('Ammonia (mmol/g/h)');
set(gca,'FontSize',FS)
