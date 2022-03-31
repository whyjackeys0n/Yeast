clear;
clc;

[DFBA_time_1,DFBA_biomass_private_1] = Yeast_MP_sucrose_function(0,1);
[DFBA_time_2,DFBA_biomass_private_2] = Yeast_MP_sucrose_function(0,2);

%% Plot
figure(1)
h1 = plot(DFBA_time_1,DFBA_biomass_private_1,'--b','linewidth',2.5);
hold on
h2 = plot(DFBA_time_2,DFBA_biomass_private_2,'--r','linewidth',2.5);
set(gca,'fontsize',24);
ylabel('Concentration [g/L]');
xlabel('Time [h]','FontSize',24);
xlim([0,DFBA_time_1(end)]);
hold on

%% load data
load CENPK21C_0611_diff_sucrose_pbpr 
CENPK21C_0611_diff_sucrose_pbpr_pop_density = CENPK21C_0611_diff_sucrose_pbpr.CENPK21C_diff_sucrose_pbpr_pop_density;
CENPK21C_0611_diff_sucrose_pbpr_time = CENPK21C_0611_diff_sucrose_pbpr.CENPK21C_diff_sucrose_pbpr_time;
CENPK21C_0611_diff_sucrose_pbpr_conc_pct = CENPK21C_0611_diff_sucrose_pbpr.CENPK21C_diff_sucrose_pbpr_conc_pct;

%% unit conversion
blank = 0.0858;
cell2gram = 4.6e-8;

%% S288C_0410_1pct_glucose_pbprch
pop_density_gL = ((CENPK21C_0611_diff_sucrose_pbpr_pop_density-blank)-0.007)./0.074.*1e6.*cell2gram;
time = CENPK21C_0611_diff_sucrose_pbpr_time;
conc_pct = CENPK21C_0611_diff_sucrose_pbpr_conc_pct;

%% Plot
pop_density_gL_mean = [];
re_exp_num = 4;
species_num = 2;
for i = 1:species_num  
    pop_density_gL_mean_unit = mean(pop_density_gL((i+1)*re_exp_num+1:(i+2)*re_exp_num,:));
    pop_density_gL_mean = [pop_density_gL_mean;pop_density_gL_mean_unit];
end
for i = 1:species_num
    plot(time,pop_density_gL_mean(i,:),'LineWidth',2);   
    hold on
end

%% Calculate interpolation
[Ttrue1,idxytrue1] = unique(DFBA_time_1);
ytrue1 = DFBA_biomass_private_1(idxytrue1);
yinterp1 = interp1(Ttrue1,ytrue1,time);
%plot(time,yinterp1,'--g','LineWidth',2);
hold on

[Ttrue2,idxytrue2] = unique(DFBA_time_2);
ytrue2 = DFBA_biomass_private_2(idxytrue2);
yinterp2 = interp1(Ttrue2,ytrue2,time);
%plot(time,yinterp2,'--y','LineWidth',2);

%% Plot settings
species = {'DFBA 1% Sucrose','DFBA 4% Sucrose','Exp Private 1% Sucrose','Exp Private 4% Sucrose'};
% species = {'DFBA 1% Sucrose','DFBA 4% Sucrose','Exp Private 1% Sucrose',...
%     'Exp Private 4% Sucrose', 'Interpolation 1% Sucrose', 'Interpolation 4% Sucrose'};
legend(species,'FontSize',24)
set(gca,'FontSize',24)
fh = figure(1);
fh.WindowState = 'maximized';
