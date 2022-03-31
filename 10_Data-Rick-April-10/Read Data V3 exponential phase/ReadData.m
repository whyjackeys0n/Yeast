clear;
clc;

S288C_1pct_glucose_pbprch = xlsread('S cerevisiae growth data from NEE.xlsx','Growth in 1% glucose','B5:CW26');
S288C_1pct_sucrose_pbprch = xlsread('S cerevisiae growth data from NEE.xlsx','Growth in 1% sucrose','B5:GX26');
S288C_diff_sucrose_pbpr = xlsread('S cerevisiae growth data from NEE.xlsx','Different sucrose concentration','B7:KO49');

CENPK21C_diff_sucrose_pb = xlsread('S cerevisiae growth data from NEE.xlsx','CEN.PK2-1C in diff. % sucrose','A4:OM16');
CENPK21C_diff_sugar_pb = xlsread('S cerevisiae growth data from NEE.xlsx','CEN.PK2-1c in diff sugars','C5:OL17');
CENPK21C_diff_sucrose_pbpr = xlsread('S cerevisiae growth data from NEE.xlsx','CEN.PK2-1c and Private unpubl.','C5:JK21');
CENPK21C_dot05pct_sucrose_pbpr = xlsread('S cerevisiae growth data from NEE.xlsx','CEN.PK2-1c and Private 0.05%S','C5:EP17');
CENPK21C_1pct_glucose_pbpr = xlsread('S cerevisiae growth data from NEE.xlsx','CEN.PK2-1c and Private 1% glu','C5:GX13');

S288C_0410_1pct_glucose_pbprch.S288C_0410_1pct_glucose_pbprch_time = S288C_1pct_glucose_pbprch(1,:);
S288C_0410_1pct_glucose_pbprch.S288C_0410_1pct_glucose_pbprch_pop_density = S288C_1pct_glucose_pbprch(2:end,:);
S288C_0410_1pct_glucose_pbprch.S288C_0410_1pct_glucose_pbprch_conc_pct = ones(size(S288C_0410_1pct_glucose_pbprch.S288C_0410_1pct_glucose_pbprch_pop_density,1),1)*1;

S288C_0410_1pct_sucrose_pbprch.S288C_0410_1pct_sucrose_pbprch_time = S288C_1pct_sucrose_pbprch(1,:);
S288C_0410_1pct_sucrose_pbprch.S288C_0410_1pct_sucrose_pbprch_pop_density = S288C_1pct_sucrose_pbprch(2:end,:);
S288C_0410_1pct_sucrose_pbprch.S288C_0410_1pct_sucrose_pbprch_conc_pct = ones(size(S288C_0410_1pct_sucrose_pbprch.S288C_0410_1pct_sucrose_pbprch_pop_density,1),1)*1;

S288C_0410_diff_sucrose_pbpr.S288C_0410_diff_sucrose_pbpr_time = S288C_diff_sucrose_pbpr(1,2:end);
S288C_0410_diff_sucrose_pbpr.S288C_0410_diff_sucrose_pbpr_pop_density = S288C_diff_sucrose_pbpr(2:end,2:end);
S288C_0410_diff_sucrose_pbpr.S288C_0410_diff_sucrose_pbpr_conc_pct = S288C_diff_sucrose_pbpr(2:end,1);

CENPK21C_0410_diff_sucrose_pb.CENPK21C_0410_diff_sucrose_pb_time = CENPK21C_diff_sucrose_pb(1,4:end);
CENPK21C_0410_diff_sucrose_pb.CENPK21C_0410_diff_sucrose_pb_pop_density = CENPK21C_diff_sucrose_pb(2:end,4:end);
CENPK21C_0410_diff_sucrose_pb.CENPK21C_0410_diff_sucrose_pb_conc_pct = CENPK21C_diff_sucrose_pb(2:end,1);
CENPK21C_0410_diff_sucrose_pb.CENPK21C_0410_diff_sucrose_pb_conc_mM = CENPK21C_diff_sucrose_pb(2:end,2);

CENPK21C_0410_diff_sugar_pb.CENPK21C_0410_diff_sugar_pb_time = CENPK21C_diff_sugar_pb(1,:);
CENPK21C_0410_diff_sugar_pb.CENPK21C_0410_diff_sugar_pb_pop_density = CENPK21C_diff_sugar_pb(2:end,:);
CENPK21C_0410_diff_sugar_pb.CENPK21C_0410_diff_sugar_pb_glucose_conc_pct = [0.5;0.5;0.5;1;1;1;0;0;0;0;0;0];
CENPK21C_0410_diff_sugar_pb.CENPK21C_0410_diff_sugar_pb_fructose_conc_pct = [0.5;0.5;0.5;0;0;0;0;0;0;0;0;0];
CENPK21C_0410_diff_sugar_pb.CENPK21C_0410_diff_sugar_pb_maltose_conc_pct = [0;0;0;0;0;0;1;1;1;0;0;0];
CENPK21C_0410_diff_sugar_pb.CENPK21C_0410_diff_sugar_pb_sucrose_conc_pct = [0;0;0;0;0;0;0;0;0;1;1;1];

CENPK21C_0410_diff_sucrose_pbpr.CENPK21C_0410_diff_sucrose_pbpr_time = CENPK21C_diff_sucrose_pbpr(1,:);
CENPK21C_0410_diff_sucrose_pbpr.CENPK21C_0410_diff_sucrose_pbpr_pop_density = CENPK21C_diff_sucrose_pbpr(2:end,:);
CENPK21C_0410_diff_sucrose_pbpr.CENPK21C_0410_diff_sucrose_pbpr_conc_pct = [1;1;1;1;4;4;4;4;1;1;1;1;4;4;4;4];

CENPK21C_0410_dot05pct_sucrose_pbpr.CENPK21C_0410_dot05pct_sucrose_pbpr_time = CENPK21C_dot05pct_sucrose_pbpr(1,:);
CENPK21C_0410_dot05pct_sucrose_pbpr.CENPK21C_0410_dot05pct_sucrose_pbpr_pop_density = CENPK21C_dot05pct_sucrose_pbpr(2:end,:);
CENPK21C_0410_dot05pct_sucrose_pbpr.CENPK21C_0410_dot05pct_sucrose_pbpr_conc_pct = ones(size(CENPK21C_0410_dot05pct_sucrose_pbpr.CENPK21C_0410_dot05pct_sucrose_pbpr_pop_density,1),1)*0.05;

CENPK21C_0410_1pct_glucose_pbpr.CENPK21C_0410_1pct_glucose_pbpr_time = CENPK21C_1pct_glucose_pbpr(1,:);
CENPK21C_0410_1pct_glucose_pbpr.CENPK21C_0410_1pct_glucose_pbpr_pop_density = CENPK21C_1pct_glucose_pbpr(2:end,:);
CENPK21C_0410_1pct_glucose_pbpr.CENPK21C_0410_1pct_glucose_pbpr_conc_pct = ones(size(CENPK21C_0410_1pct_glucose_pbpr.CENPK21C_0410_1pct_glucose_pbpr_pop_density,1),1)*1;

save S288C_0410_1pct_glucose_pbprch S288C_0410_1pct_glucose_pbprch
save S288C_0410_1pct_sucrose_pbprch S288C_0410_1pct_sucrose_pbprch
save S288C_0410_diff_sucrose_pbpr S288C_0410_diff_sucrose_pbpr
save CENPK21C_0410_diff_sucrose_pb CENPK21C_0410_diff_sucrose_pb
save CENPK21C_0410_diff_sugar_pb CENPK21C_0410_diff_sugar_pb
save CENPK21C_0410_diff_sucrose_pbpr CENPK21C_0410_diff_sucrose_pbpr
save CENPK21C_0410_dot05pct_sucrose_pbpr CENPK21C_0410_dot05pct_sucrose_pbpr
save CENPK21C_0410_1pct_glucose_pbpr CENPK21C_0410_1pct_glucose_pbpr












