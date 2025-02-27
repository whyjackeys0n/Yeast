%% clear environment
clear;
clc;

% Sheet 1 - S288C_0410_1pct_glucose_pbprch
% Sheet 2 - S288C_0410_1pct_sucrose_pbprch
% Sheet 3 - S288C_0410_diff_sucrose_pbpr

% Sheet 4 - CENPK21C_0410_diff_sucrose_pb
% Sheet 5 - CENPK21C_0410_diff_sugar_pb
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
pop_density_gL_mean = [];
re_exp_num = 7;
species_num = 3;
color_list = ['b','r','y','m'];

figure(1)
subplot(1,2,1)
for i = 1:species_num 
    for j = re_exp_num*(i-1)+1:re_exp_num*i
        plot(time,log(pop_density_gL(j,:)),'LineWidth',2,'Color',color_list(i));
        hold on
    end
end 
set(gca,'FontSize',24)

for i = 1:species_num   
    pop_density_gL_mean_unit = zeros(1,length(pop_density_gL));
    for j = re_exp_num*(i-1)+1:re_exp_num*(i-1)+re_exp_num
        pop_density_gL_mean_unit = pop_density_gL_mean_unit + pop_density_gL(j,:);
    end
    pop_density_gL_mean = [pop_density_gL_mean;pop_density_gL_mean_unit./re_exp_num];
end
subplot(1,2,2)
for i = 1:species_num
    plot(time,log(pop_density_gL_mean(i,:)),'LineWidth',2);    
    hold on
end
legend(species,'FontSize',24)
set(gca,'FontSize',24)

%% fitting
[~,idx_start] = max(time>=3,[],2);
[~,idx_end] = max(time>=13,[],2);

% time scale
time0 = time(idx_start:idx_end)-time(idx_start);

figure(2)
for i = 1:species_num
    subplot(1,species_num,i)
    x0 = pop_density_gL_mean(i,idx_start);
    lnxtx0 = log(pop_density_gL_mean(i,idx_start:idx_end)./x0);
    p{i} = regress(lnxtx0',time0');
    y = time0.*p{i};
    plot(time0,y,time0,lnxtx0,'LineWidth',2);
    set(gca,'FontSize',24)
    disp([species{i},': ',num2str(p{i})]);
end





