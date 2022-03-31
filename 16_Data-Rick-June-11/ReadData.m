clear;
clc;

CENPK21C_diff_sucrose_pbpr = xlsread...
    ('Private metaboliser and public metaboliser CEN.PK2-1C 1 and 4% sucrose OD.xlsx',...
    'Private and Public 1% 4% Sucros','C4:JK20');

CENPK21C_0611_diff_sucrose_pbpr.CENPK21C_diff_sucrose_pbpr_time = CENPK21C_diff_sucrose_pbpr(1,:);
CENPK21C_0611_diff_sucrose_pbpr.CENPK21C_diff_sucrose_pbpr_pop_density = CENPK21C_diff_sucrose_pbpr(2:end,:);
CENPK21C_0611_diff_sucrose_pbpr.CENPK21C_diff_sucrose_pbpr_conc_pct = [1;1;1;1;4;4;4;4;1;1;1;1;4;4;4;4];

save CENPK21C_0611_diff_sucrose_pbpr CENPK21C_0611_diff_sucrose_pbpr













