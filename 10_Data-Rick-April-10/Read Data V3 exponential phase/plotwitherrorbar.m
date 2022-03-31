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
load S288C_0410_1pct_sucrose_pbprch
load S288C_0410_diff_sucrose_pbpr

%% unit convert from cell/microliter to gram/liter
unit = 1e6 * 4.6e-11;

%% S288C_0410_1pct_glucose_pbprch
pop_density_gL = S288C_0410_1pct_glucose_pbprch.S288C_0410_1pct_glucose_pbprch_pop_density.*unit;
time = S288C_0410_1pct_glucose_pbprch.S288C_0410_1pct_glucose_pbprch_time;
conc_pct = S288C_0410_1pct_glucose_pbprch.S288C_0410_1pct_glucose_pbprch_conc_pct;
species = {'Private','Cheat','Public'};

%% figure
re_exp_num = 7;
species_num = 3;

figure(1)
for i = 1:species_num
    subplot(1,species_num,i)
    [mu,sigma] = normfit(pop_density_gL(re_exp_num*(i-1)+1:re_exp_num*i,:),0.05);
    plot(time,mu,'LineWidth',2,'Color','b');
    hold on
    plot(time,mu+sigma,'.','LineWidth',2,'Color','r');
    hold on
    plot(time,mu-sigma,'.','LineWidth',2,'Color','r');
    legend(species{i},'FontSize',24)
    set(gca,'Fontsize',24)
end






