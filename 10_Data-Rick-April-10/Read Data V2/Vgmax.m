%% clear environment
clear;
clc;

% Sheet 1 - S288C_0410_1pct_glucose_pbprch
% Sheet 2 - S288C_0410_1pct_sucrose_pbprch
% Sheet 3 - S288C_0410_diff_sucrose_pbpr

% Sheet 4 - CENPK21C_0410_diff_sucrose_pb
% Sheet 5 - CENPK21C_0410_1pct_sugar_pb
% Sheet 6 - CENPK21C_0410_diff_sucrose_pbpr
% Sheet 7 - CENPK21C_0410_dot05pct_sucrose_pbpr
% Sheet 8 - CENPK21C_0410_1pct_glucose_pbpr

%% load data
load S288C_0410_1pct_glucose_pbprch
S288C_0410_1pct_glucose_pbprch_pop_density = S288C_0410_1pct_glucose_pbprch.S288C_0410_1pct_glucose_pbprch_pop_density;
S288C_0410_1pct_glucose_pbprch_time = S288C_0410_1pct_glucose_pbprch.S288C_0410_1pct_glucose_pbprch_time;
S288C_0410_1pct_glucose_pbprch_conc_pct = S288C_0410_1pct_glucose_pbprch.S288C_0410_1pct_glucose_pbprch_conc_pct;

%% unit convert from cell/microliter to gram/liter
unit = 1e6 * 4.6e-11;

%% S288C_0410_1pct_glucose_pbprch
pop_density_gL = S288C_0410_1pct_glucose_pbprch_pop_density.*unit;
time = S288C_0410_1pct_glucose_pbprch_time;
conc_pct = S288C_0410_1pct_glucose_pbprch_conc_pct;
species = {'Private','Cheat','Public'};

%% figure 1
figure(1)
subplot(1,2,1)
for i = 1:length(conc_pct)
    plot(time,log(pop_density_gL(i,:)),'LineWidth',2);
    hold on
end

pop_density_gL_mean = [];
for i = 1:3   
    pop_density_gL_mean_unit = zeros(1,length(pop_density_gL));
    for j = (i-1)*7+1:(i-1)*7+7
        pop_density_gL_mean_unit = pop_density_gL_mean_unit + pop_density_gL(j,:);
    end
    pop_density_gL_mean = [pop_density_gL_mean;pop_density_gL_mean_unit./7];
end
subplot(1,2,2)
for i = 1:3
    plot(time,log(pop_density_gL_mean(i,:)),'LineWidth',2);    
    hold on
end
legend(species,'FontSize',20)
set(gca,'FontSize',20)

%% fitting
[~,idx_start] = max(time>=3,[],2);
[~,idx_end] = max(time>=13,[],2);

% time scale
time0 = time(idx_start:idx_end)-time(idx_start);

figure(2)
subplot(1,3,1)
x0 = pop_density_gL_mean(1,idx_start);
lnxtx0 = log(pop_density_gL_mean(1,idx_start:idx_end)./x0);
% p{1} = polyfit(time0,lnxtx0,1);
% y = polyval(p{1},time0);
p{1} = regress(lnxtx0',time0');
y = time0.*p{1};
plot(time0,y,time0,lnxtx0,'LineWidth',2);

subplot(1,3,2)
x0 = pop_density_gL_mean(2,idx_start);
lnxtx0 = log(pop_density_gL_mean(2,idx_start:idx_end)./x0);
% p{2} = polyfit(time0,lnxtx0,1);
% y = polyval(p{2},time0);
p{2} = regress(lnxtx0',time0');
y = time0.*p{2};
plot(time0,y,time0,lnxtx0,'LineWidth',2);

subplot(1,3,3)
x0 = pop_density_gL_mean(3,idx_start);
lnxtx0 = log(pop_density_gL_mean(3,idx_start:idx_end)./x0);
% p{3} = polyfit(time0,lnxtx0,1);
% y = polyval(p{3},time0);
p{3} = regress(lnxtx0',time0');
y = time0.*p{3};
plot(time0,y,time0,lnxtx0,'LineWidth',2);

disp(['private growth rate: ',num2str(p{1})]);
disp(['cheat growth rate: ',num2str(p{2})]);
disp(['public growth rate: ',num2str(p{3})]);





