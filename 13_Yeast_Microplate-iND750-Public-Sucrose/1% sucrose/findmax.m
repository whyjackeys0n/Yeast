%% Program for find the exponential phase ending point
clear;
clc;

%% Load data
load loop52018510050400

%% Set the expect range
xlb = 19.8177;
xub = 20.8427;
ylb = 3.29676;
yub = 3.42996;
figure(100)
plot(Tn_unit,Xn_unit,'.');
hold on
rectangle('position',[xlb,ylb,xub-xlb,yub-ylb],'LineWidth',2)

%% load data
load CENPK21C_0410_diff_sucrose_pbpr 
CENPK21C_0410_diff_sucrose_pbpr_pop_density = CENPK21C_0410_diff_sucrose_pbpr.CENPK21C_0410_diff_sucrose_pbpr_pop_density;
CENPK21C_0410_diff_sucrose_pbpr_time = CENPK21C_0410_diff_sucrose_pbpr.CENPK21C_0410_diff_sucrose_pbpr_time;
CENPK21C_0410_diff_sucrose_pbpr_conc_pct = CENPK21C_0410_diff_sucrose_pbpr.CENPK21C_0410_diff_sucrose_pbpr_conc_pct;

%% unit convert from cell/microliter to gram/liter
unit = 1e6 * 4.6e-11;

%% S288C_0410_1pct_glucose_pbprch
pop_density_gL = CENPK21C_0410_diff_sucrose_pbpr_pop_density.*unit;
time = CENPK21C_0410_diff_sucrose_pbpr_time;
conc_pct = CENPK21C_0410_diff_sucrose_pbpr_conc_pct;
species = {'Sim','Exp Public','Exp Private'};

%% figure 1
pop_density_gL_mean = [];
re_exp_num = 4;
species_num = 2;

for i = 1:species_num   
    pop_density_gL_mean_unit = zeros(1,length(pop_density_gL));
    for j = re_exp_num*(i-1)+1:re_exp_num*(i-1)+re_exp_num
        pop_density_gL_mean_unit = pop_density_gL_mean_unit + pop_density_gL(j,:);
    end
    pop_density_gL_mean = [pop_density_gL_mean;pop_density_gL_mean_unit./re_exp_num];
end
for i = 1:species_num
    plot(time,pop_density_gL_mean(i,:),'LineWidth',2);   
    hold on
end
set(gca,'FontSize',24)

%% Find parameters
idxTn = find(Tn_unit>xlb&Tn_unit<xub);
idxXn = find(Xn_unit>ylb&Xn_unit<yub);
idxTnXn = intersect(idxTn,idxXn);
figure(200)
plot(Tn_unit(idxTnXn),Xn_unit(idxTnXn),'o','LineWidth',5);
set(gca,'FontSize',24)

%% Display parameter sets
steploop = 0;
stepres = 1;
res_unit = [];
res = [];
for i = 5:1:20
    for j = 1:1:8
        for k = 5:5:100
            for l = 50:50:400
                steploop = steploop + 1;
                for m = 1:length(idxTnXn)
                    if steploop == idxTnXn(m)
                        res_unit = [i,j,k,l];
                        res = [res;res_unit];
                        fitMSE(stepres) = Yeast_plot_function(stepres,i,j,k,l);
                        stepres = stepres + 1;
                    end
                end
            end
        end
    end
end

for n = 1:stepres-1
    disp(['Vgmax',num2str(res(n,1)),' Vomax',num2str(res(n,2)),' KLa',num2str(res(n,3)),' Xi',num2str(res(n,4))]);
end
