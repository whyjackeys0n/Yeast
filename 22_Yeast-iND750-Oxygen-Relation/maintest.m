clear;
clc;

grpb_res = [];
for i = 1:0.1:10
    grpb = testoxygen(i,'iND750');
    grpb_res = [grpb_res,grpb];
end

grpr_res = [];
for i = 1:0.1:10
    grpr = testoxygen(i,'iND750_Private');
    grpr_res = [grpr_res,grpr];
end

grch_res = [];
for i = 1:0.1:10
    grch = testoxygen(i,'iND750_Cheat');
    grch_res = [grch_res,grch];
end

plot(1:0.1:10,grpb_res,1:0.1:10,grpr_res);

