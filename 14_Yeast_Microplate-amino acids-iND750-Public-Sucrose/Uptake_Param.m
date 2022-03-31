%% Script and Source of the uptake parameters
% J. Wang, 3/22/2021
clear;
clc;

%% Glucose
% Jared L. Hjersted, Michael A. Henson, Radhakrishnan Mahadevan. 2007.
% Genome-Scale Analysis of Saccharomyces cerevisiae Metabolism and Ethanol 
% Production in Fed-Batch Culture. Biotechnology and Bioengineering. 97:
% 1190-1204.

% B. Sonnleitnert and O. Kappeli. 1986. Growth of Saccharomyces cerevisiae 
% Is Controlled by Its Limited Respiratory Capacity: Formulation and 
% Verification of a Hypothesis. Biotechnology and Bioengineering. 28:
% 927-937.

qsmax = 3.5;    % maximal specific glucose uptake rate (g·g-1·h-1)
Ks = 0.3;       % saturation parameter for glucose uptake (g·L-1)

VGmax_ref = qsmax/180*1000;    
KmG_ref = Ks/180*1000;

VGmax = VGmax_ref;     % maximum glucose uptake rate (mmol/gDW/h)
KmG = KmG_ref;      % glucose uptake saturation constants (gDW/L)

clear qsmax Ks

%% Fructose
% Carole Guillaume, Pierre Delobel, Jean-Marie Sablayrolles, Bruno Blondin.
% 2007. Molecular Basis of Fructose Utilization by the Wine Yeast 
% Saccharomyces cerevisiae: a Mutated HXT3 Allele Enhances Fructose 
% Fermentation. Applied and Environmental Microbiology. 73: 2432-2439.

Vmax_glucose = 210;     % Hxt3-V5 maximum glucose uptake rate (nmol·min-1·mg[dry wt]-1)
Vmax_fructose = 273;    % Hxt3-V5 maximum fructose uptake rate (nmol·min-1·mg[dry wt]-1)
Km_glucose = 2.1;       % high-affinity glucose transport kinetics (mM)
Km_fructose = 4.6;      % high-affinity fructose transport kinetics (mM)

VFmax_ref = VGmax*Vmax_fructose/Vmax_glucose;
KmF_ref = KmG*Km_fructose/Km_glucose;

VFmax = 26;
KmF = 1;

clear Vmax_glucose Vmax_fructose Km_glucose Km_fructose

%% Sucrose
% Ramon Wahl1, Kathrin Wippel, Sarah Goos, Jorg Kamper, Norbert Sauer.
% 2010. A Novel High-Affinity Sucrose Transporter Is Required for Virulence
% of the Plant Pathogen Ustilago maydis. PLoS Biology. 8: e1000303.

% Frans M. Klis, Chris G. de Koster, Stanley Brul. 2014. Cell Wall-Related 
% Bionumbers and Bioestimates of Saccharomyces cerevisiae and Candida 
% albicans. Eukaryotic Cell. 13: 2-9.

Vmax_sucrose =0.77;% Maximum fructose uptake rate (nmoles·min-1·mg FW-1)
Km_sucrose = 26;    % Michaelis-Menten kinetics of sucrose uptake rates (μM)
DW_frac = 0.34;     % dry weight fraction of wet weight

VSmax_ref = Vmax_sucrose*1e-6*60*1000/DW_frac;
KmS_ref = Km_sucrose/1000;

VSmax = 0.136;
KmS = 0.026;

clear Vmax_sucrose Km_sucrose DW_frac



