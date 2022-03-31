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

figure(1)
subplot(1,2,1)
for i = 1:length(conc_pct)
    plot(time,pop_density_gL(i,:));
    hold on
end

pop_density_mean_gL = [];
for i = 1:3   
    pop_density_gL_mean_unit = zeros(1,length(pop_density_gL));
    for j = (i-1)*7+1:(i-1)*7+7
        pop_density_gL_mean_unit = pop_density_gL_mean_unit + pop_density_gL(j,:);
    end
    pop_density_mean_gL = [pop_density_mean_gL;pop_density_gL_mean_unit./7];
end
subplot(1,2,2)
for i = 1:3  
    plot(time,pop_density_mean_gL(i,:));
    hold on
end
legend(species,'FontSize',20)

%% S288C_0410_1pct_sucrose_pbprch
clear pop_density_gL time conc_pct
pop_density_gL = S288C_0410_1pct_sucrose_pbprch.S288C_0410_1pct_sucrose_pbprch_pop_density.*unit;
time = S288C_0410_1pct_sucrose_pbprch.S288C_0410_1pct_sucrose_pbprch_time;
conc_pct = S288C_0410_1pct_sucrose_pbprch.S288C_0410_1pct_sucrose_pbprch_conc_pct;
species = {'Private','Cheat','Public'};

figure(2)
subplot(1,2,1)
for i = 1:length(conc_pct)
    plot(time,pop_density_gL(i,:));
    hold on
end

pop_density_mean_gL = [];
for i = 1:3   
    pop_density_gL_mean_unit = zeros(1,length(pop_density_gL));
    for j = (i-1)*7+1:(i-1)*7+7
        pop_density_gL_mean_unit = pop_density_gL_mean_unit + pop_density_gL(j,:);
    end
    pop_density_mean_gL = [pop_density_mean_gL;pop_density_gL_mean_unit./7];
end
subplot(1,2,2)
for i = 1:3  
    plot(time,pop_density_mean_gL(i,:));   
    hold on
end
legend(species,'FontSize',20)

%% S288C_0410_diff_sucrose_pbpr
clear pop_density_gL time conc_pct
pop_density_gL = S288C_0410_diff_sucrose_pbpr.S288C_0410_diff_sucrose_pbpr_pop_density.*unit;
time = S288C_0410_diff_sucrose_pbpr.S288C_0410_diff_sucrose_pbpr_time;
conc_pct = S288C_0410_diff_sucrose_pbpr.S288C_0410_diff_sucrose_pbpr_conc_pct;
species = {'0.0625','0.125','0.25','0.5','1','2','4'};

figure(3)
subplot(1,2,1)
for i = 1:21
    plot(time,pop_density_gL(i,:));
    hold on
end

pop_density_mean_gL = [];
for i = 1:7   
    pop_density_gL_mean_unit = zeros(1,length(pop_density_gL));
    for j = (i-1)*3+1:(i-1)*3+3
        pop_density_gL_mean_unit = pop_density_gL_mean_unit + pop_density_gL(j,:);
    end
    pop_density_mean_gL = [pop_density_mean_gL;pop_density_gL_mean_unit./7];
end
subplot(1,2,2)
for i = 1:7  
    plot(time,pop_density_mean_gL(i,:));
    hold on
end
legend(species,'FontSize',20)
