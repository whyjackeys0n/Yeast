clear;
clc;

[DFBA_time_1,DFBA_biomass_private_1] = Yeast_MP_sucrose_function(0.0625,1);
[DFBA_time_2,DFBA_biomass_private_2] = Yeast_MP_sucrose_function(0.25,2);

%% Plot
figure(1)
h1 = plot(DFBA_time_1,DFBA_biomass_private_1,'linewidth',2.5);
hold on
h2 = plot(DFBA_time_2,DFBA_biomass_private_2,'linewidth',2.5);
set(gca,'fontsize',30);
ylabel('Concentration [g/L]');
xlabel('Time [h]','FontSize',30);
% xlim([0,DFBA_time_1(end)]);
hold on

%% load data
load CENPK21C_0410_diff_sucrose_pb 
CENPK21C_0410_diff_sucrose_pb_pop_density = CENPK21C_0410_diff_sucrose_pb.CENPK21C_0410_diff_sucrose_pb_pop_density;
CENPK21C_0410_diff_sucrose_pb_time = CENPK21C_0410_diff_sucrose_pb.CENPK21C_0410_diff_sucrose_pb_time;
CENPK21C_0410_diff_sucrose_pb_conc_pct = CENPK21C_0410_diff_sucrose_pb.CENPK21C_0410_diff_sucrose_pb_conc_pct;

%% S288C_0410_1pct_glucose_pbprch
pop_density_gL = CENPK21C_0410_diff_sucrose_pb_pop_density*1e6*4.6e-11/5.2;
time = CENPK21C_0410_diff_sucrose_pb_time;
conc_pct = CENPK21C_0410_diff_sucrose_pb_conc_pct;

%% Plot
pop_density_gL_mean = [];
re_exp_num = 3;
species_num = 2;
species_start_num = 3;
for i = species_start_num:species_start_num+species_num-1  
    pop_density_gL_mean_unit = mean(pop_density_gL((i-1)*re_exp_num+1:i*re_exp_num,:));
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
species = {'DFBA 0.0625% Sucrose','DFBA 0.25% Sucrose'};
% species = {'DFBA 1% Sucrose','DFBA 4% Sucrose','Exp Private 1% Sucrose',...
%     'Exp Private 4% Sucrose', 'Interpolation 1% Sucrose', 'Interpolation 4% Sucrose'};
legend(species,'FontSize',30)
set(gca,'FontSize',30)
fh = figure(1);
fh.WindowState = 'maximized';
