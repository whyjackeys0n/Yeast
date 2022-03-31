clear;
clc;

load('iMM904.mat');

%% Private Model
model_private = iMM904;

idx_scr_hydro = find(contains(model_private.rxnNames,'Sucrose hydrolyzing enxyme  extracellular'));
idx_scr_trans = find(contains(model_private.rxnNames,'Sucrose transport in via proton symport'));

idx_ATP = find(strcmp(model_private.mets,'atp_c'));
idx_ADP = find(strcmp(model_private.mets,'adp_c'));
idx_Pi = find(strcmp(model_private.mets,'pi_c'));
% idx_H2O = find(strcmp(model.mets,'h2o_c'));
% idx_H = find(strcmp(model.mets,'h_c'));

model_private.S([idx_ATP],idx_scr_trans) = -1;
model_private.S([idx_ADP,idx_Pi],idx_scr_trans) = 1;

model_private.S(:,idx_scr_hydro) = 0;
model_private.S(620,idx_scr_hydro) = -1;
model_private.S(1015,idx_scr_hydro) = -1;
model_private.S(592,idx_scr_hydro) = 1;
model_private.S(506,idx_scr_hydro) = 1;

model_private.rxnNames{idx_scr_hydro} = 'Sucrose hydrolyzing enxyme  cytosolic';

iMM904_Private = model_private;
save iMM904_Private iMM904_Private

%% Cheat Model
model_cheat = iMM904;
model_cheat.S(:,idx_scr_hydro) = 0;
model_cheat.S(:,idx_scr_trans) = 0;

iMM904_Cheat = model_cheat;
save iMM904_Cheat iMM904_Cheat


