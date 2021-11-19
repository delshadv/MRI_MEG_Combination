% Main Script:
% This is the "main" script to reproduce results of "Late Combination shows MEG adds 
% to MRI in classifying MCI versus controls" 
% published in ""
% This script carries out permutation test for MCI vs Control classification
% problem using combination of MEG - MRI  with Multi-Kernel Learning.

% You might need to use Parallel Computing Toolbox to be able to run MKL functions
% in parallel for many repetitions to get more accurate results.
%
% by Rik Henson and Delshad Vaghari, 2020

% TO DO:    re-run on single broadband power - Done

%% Define Path to needed functions

% Assumed you are currently in the directory including BioFIND data,
% OSL and MKL directories as described in readme.md
bwd = pwd; % needs to be in Github directory
addpath(fullfile(bwd,'MKL'));
addpath(fullfile(bwd,'MKL','supplementary'));
addpath(fullfile(bwd,'MKL','derived'));

cd MKL

%% Analysis 1 (Comparison of MRI vs MEG vs COFs through all possible combination methods)
% Import data and define input cell array
MRI = csvread('derived/ROIdata.csv');
MEG = csvread('derived/GRDCOVGAMMA1.csv'); % Since lgamma does best numerically
labels = csvread('derived/labels.csv');
COFS = csvread('derived/COFS.csv');

COFS_cell = {COFS(:,1),COFS(:,2),COFS(:,3),COFS(:,4),COFS(:,4),COFS(:,6),COFS(:,7),COFS(:,8)};

V = {{COFS(:,1)},COFS_cell,{MRI},{MEG},[COFS_cell(:);MRI]',[COFS_cell(:);{MEG}]',{MRI,MEG},...
    [COFS_cell(:);MRI;MEG]'};

% Tune Hyper-parameters
rng('default') % For reproducibility
[bLcomI,bLcomL,~,~] = mkl_ens_nestedcv(V,labels,'Nfold_o',5,'Nfold_i',10,'Nrun',10);

LI = mean(mean(bLcomI,3)); % Lambda for Intermediate
% 0.5704   0.0100   0.1606    0.4000    0.1792    0.0914    0.1450    0.1484
LL = mean(mean(bLcomL,3)); % Lambda for Late
% 1.0000   0.6036    1.0000    1.0000    0.5300    0.6956    0.7416   0.6734

S = struct;
S.x1 = LI(1); % 0.5704 hyperparameter for COF only
S.x1dim = 1;  % Dimension of COF only
S.x2 = LI(3); % 0.1606 hyperparameter for MRI only
S.x2dim = size(MRI,2); % Dimension of MRI only
S.x3 = LI(4); % 0.4000; hyperparameter for MRI only
S.x3dim = size(MEG,2); % Dimension of MEG only

% Classification step
V = {COFS_cell,{MRI},{MEG},[COFS_cell(:);MRI]',[COFS_cell(:);{MEG}]',{MRI,MEG},...
    [COFS_cell(:);MRI;MEG]'};

LI(1) = [];
LL(1) = [];

rng('default') % For reproducibility
[acc1,acc2] = mkl_ens_n(V,labels,'Nfold',5,'Nrun',1000,'feat_norm',1,'ens',1,...
    'machine','kernel','HyperI',LI,'HyperL',LL,...
    'Hyper_sizeS',S);

clear V
V = {{COFS},{MRI},{MEG},{[COFS MRI]},{[COFS MEG]},{[MRI MEG]},{[COFS MRI MEG]}};

% Tune Hyper-parameters
rng('default') % For reproducibility
[bLcomI,~,~,~] = mkl_ens_nestedcv(V,labels,'Nfold_o',5,'Nfold_i',10,'Nrun',10,'ens',0);

LI = mean(mean(bLcomI,3)); % Lambda for Intermediate
% 0.4294    0.1756    0.4442    0.1466    0.4000    0.3046    0.3472

rng('default') % For reproducibility
acc3 = mkl_ens_n(V,labels,'Nfold',5,'Nrun',1000,'feat_norm',1,'ens',0,...
    'machine','kernel','HyperI',LI);

save FSCofsMRIMEG_GrdCov acc1 acc2 acc3

% Plot resluts (Classification accuracy and Pos-hoc comparison)

% define titles in same order of the input cell V
titles = {'COF','MRI','MEG','MRI,COF','MEG,COF'};
pos_titles = {'MRI>COF','MEG>COF','MEG>MRI',...
    'MRI,COF>COF','MEG,COF>COF'};
%  define contrasts
c = [-1 1 0  0 0 ;
    -1 0 1   0 0 ;
    0 -1 1   0 0 ;
    
    -1 0 0   1 0 ;
    -1 0 0   0 1 ];

f1 = plot_results(0,titles,acc2(:,1:5,:),pos_titles,c); % main figure 1
sgtitle('Late Combination')
sgt.FontSize = 20;
% Save Figures
saveas(f1,'Fig1.png')

f2 = plot_results(0,titles,acc3(:,1:5,:),pos_titles,c); % Supplementary figure 3
sgtitle('Early Combination')
sgt.FontSize = 20;
% Save Figures
saveas(f2,'Supp Fig3.png')

f3 = plot_results(0,titles,acc1(:,1:5,:),pos_titles,c); % Supplementary figure 4
sgtitle('Intermediate Combination')
sgt.FontSize = 20;

% Save Figures
saveas(f3,'Supp Fig4.png')

accuracy = [];
accuracy = cat(2,acc3(:,6,:),acc1(:,6,:));
accuracy = cat(2,accuracy,acc2(:,[6 7 2 4],:));

% Plot resluts (Classification accuracy and Pos-hoc comparison)

% define titles in same order of the input cell V
titles = {'MRI,MEG (Early Comb.)','MRI,MEG (Intermediate Comb.)',...
    'MRI,MEG (Late Comb.)','MRI,MEG,COF (Late Comb.)'};%,'MRI,COF (Late)','MRI,MEG,COF (Late)'

pos_titles = {'MRI,MEG (Intermediate Comb.) > MRI,MEG (Early Comb.)',...
    'MRI,MEG (Late Comb.) > MRI,MEG (Intermediate Comb.)',...
    'MRI,MEG (Late Comb.) > MRI','MRI,MEG,COF (Late Comb.) > MRI,COF (Late Comb.)'};

%  define contrasts
c = [-1 1 0  0 0 0;
    0 -1 1  0 0 0;
    0 0 1  0 -1 0;
    0 0 0  1 0 -1];

f4 = plot_results(0,titles,accuracy(:,1:end-2,:),pos_titles,c,accuracy);
sgtitle('Early vs Intermediate vs Late Combination')
sgt.FontSize = 20;
saveas(f4,'Fig2.png')

%% Analysis  (Other MEG - MRI combination), MEG: GRD COV
% Import data and define input cell array
clear
tbl = table();
tbl.bands = {'Delta(2-4Hz)';'Theta(4-8Hz)';'Alpha(8-12Hz)';'Beta(12-30Hz)';'Low-Gamma(30-48)';'High-Gamma(52-86)'};

MRI = csvread('derived/ROIdata.csv');
labels = csvread('derived/labels.csv');

V = {{csvread('GRDCOVDELTA.csv')},{csvread('GRDCOVTHETA.csv')},{csvread('GRDCOVALPHA.csv')}...
    ,{csvread('GRDCOVBETA.csv')},{csvread('GRDCOVGAMMA1.csv')},{csvread('GRDCOVGAMMA2.csv')},...
    {MRI,csvread('GRDCOVDELTA.csv')},{MRI,csvread('GRDCOVTHETA.csv')},{MRI,csvread('GRDCOVALPHA.csv')}...
    ,{MRI,csvread('GRDCOVBETA.csv')},{MRI,csvread('GRDCOVGAMMA1.csv')},{MRI,csvread('GRDCOVGAMMA2.csv')}};

% Tune Hyper-parameters
rng('default') % For reproducibility
[bLcomI,bLcomL,~,~] = mkl_ens_nestedcv(V,labels,'Nfold_o',5,'Nfold_i',10,'Nrun',10);

LI = mean(mean(bLcomI,3)); % Lambda for Intermediate
% 0.1624    0.2892    0.2862    0.3052    0.4490    0.4360    0.2888
% 0.2188    0.1652    0.1332 0.1406    0.1344
LL = mean(mean(bLcomL,3)); % Lambda for Late
% 1.0000    1.0000    1.0000    1.0000    1.0000    1.0000    0.3852
% 0.4520    0.6338    0.7130    0.7094    0.5164

V = {{MRI},{MRI,csvread('GRDCOVDELTA.csv')},{MRI,csvread('GRDCOVTHETA.csv')},{MRI,csvread('GRDCOVALPHA.csv')}...
    ,{MRI,csvread('GRDCOVBETA.csv')},{MRI,csvread('GRDCOVGAMMA1.csv')},{MRI,csvread('GRDCOVGAMMA2.csv')}};

S = struct;
S.x1 = 0.1606; % hyperparameter for MRI only (taken from the first analysis)
S.x2 = LI(1);
S.x3 = LI(2);
S.x4 = LI(3);
S.x5 = LI(4);
S.x6 = LI(5);
S.x7 = LI(6);
S.x1dim = size(MRI,2); % Dimension of MRI only
S.x2dim = size(V{1,2}{1,2},2);  % Dimension of MEG only (delta)
S.x3dim = size(V{1,2}{1,2},2);  % Dimension of MEG only (theta)
S.x4dim = size(V{1,2}{1,2},2);  % Dimension of MEG only (alpha)
S.x5dim = size(V{1,2}{1,2},2);  % Dimension of MEG only (beta)
S.x6dim = size(V{1,2}{1,2},2);  % Dimension of MEG only (gamma1)
S.x7dim = size(V{1,2}{1,2},2);  % Dimension of MEG only (gamma2)

% Classification step
LI = [0.1606 LI(7) LI(8) LI(9) LI(10) LI(11) LI(12)];
LL = [1.0000 LL(7) LL(8) LL(9) LL(10) LL(11) LL(12)];

rng('default') % For reproducibility
[~,acc2] = mkl_ens_n(V,labels,'Nfold',5,'Nrun',1000,'feat_norm',1,'ens',1,...
    'machine','kernel','HyperI',LI,'HyperL',LL,...
    'Hyper_sizeS',S);

save MRIMEG_GrdCov acc2
tbl.COV_of_GRD = [mean(mean(acc2(:,2:end,:),3))' std(mean(acc2(:,2:end,:),3))'];

% Plot resluts (Classification accuracy and Pos-hoc comparison)

% define titles in same order of the input cell V
titles = {'MRI','MRI,MEG(GRD COV Delta)','MRI,MEG(GRD COV Theta)','MRI,MEG(GRD COV Alpha)',...
    'MRI,MEG(GRD COV Beta)','MRI,MEG(GRD COV lGamma)','MRI,MEG(GRD COV hGamma)'};

pos_titles = {'MRI,MEG(GRD COV Delta)>MRI','MRI,MEG(GRD COV Theta)>MRI','MRI,MEG(GRD COV Alpha)>MRI',...
    'MRI,MEG(GRD COV Beta)>MRI','MRI,MEG(GRD COV lGamma)>MRI','MRI,MEG(GRD COV hGamma)>MRI'};

%  define contrasts
c = [-1 1 0  0 0 0  0;
    -1 0 1  0 0 0  0;
    -1 0 0  1 0 0  0;
    
    -1 0 0  0 1 0  0;
    -1 0 0  0 0 1  0;
    -1 0 0  0 0 0  1];

f5 = plot_results(1,titles,acc2,pos_titles,c);
sgtitle('MEG (GRD COV) MRI Combinations')
sgt.FontSize = 20;

%% Analysis  (Other MEG - MRI combination), MEG: GRD VAR
% Import data and define input cell array
clear
MRI = csvread('derived/ROIdata.csv');
labels = csvread('derived/labels.csv');

V = {{csvread('GRDVARDELTA.csv')},{csvread('GRDVARTHETA.csv')},{csvread('GRDVARALPHA.csv')},...
    {csvread('GRDVARBETA.csv')},{csvread('GRDVARGAMMA1.csv')},{csvread('GRDVARGAMMA2.csv')},...
    {MRI,csvread('GRDVARDELTA.csv')},{MRI,csvread('GRDVARTHETA.csv')},{MRI,csvread('GRDVARALPHA.csv')}...
    ,{MRI,csvread('GRDVARBETA.csv')},{MRI,csvread('GRDVARGAMMA1.csv')},{MRI,csvread('GRDVARGAMMA2.csv')}};

% Tune Hyper-parameters
rng('default') % For reproducibility
[bLcomI,bLcomL,~,~] = mkl_ens_nestedcv(V,labels,'Nfold_o',5,'Nfold_i',10,'Nrun',10);

LI = mean(mean(bLcomI,3)); % Lambda for Intermediate
% 0.2822    0.3442    0.3168    0.3076    0.2172    0.2968    0.3938
% 0.4044    0.2024    0.1952    0.1400    0.1476
LL = mean(mean(bLcomL,3)); % Lambda for Late
% 1.0000    1.0000    1.0000    1.0000    1.0000    1.0000    0.3642
% 0.4328    0.6306    0.6250    0.5244    0.5104

V = {{MRI},{MRI,csvread('GRDVARDELTA.csv')},{MRI,csvread('GRDVARTHETA.csv')},{MRI,csvread('GRDVARALPHA.csv')}...
    ,{MRI,csvread('GRDVARBETA.csv')},{MRI,csvread('GRDVARGAMMA1.csv')},{MRI,csvread('GRDVARGAMMA2.csv')}};

S = struct;
S.x1 = 0.1606; % hyperparameter for MRI only (taken from the first analysis)
S.x2 = LI(1);
S.x3 = LI(2);
S.x4 = LI(3);
S.x5 = LI(4);
S.x6 = LI(5);
S.x7 = LI(6);
S.x1dim = size(MRI,2); % Dimension of MRI only
S.x2dim = size(V{1,2}{1,2},2);  % Dimension of MEG only (delta)
S.x3dim = size(V{1,2}{1,2},2);  % Dimension of MEG only (theta)
S.x4dim = size(V{1,2}{1,2},2);  % Dimension of MEG only (alpha)
S.x5dim = size(V{1,2}{1,2},2);  % Dimension of MEG only (beta)
S.x6dim = size(V{1,2}{1,2},2);  % Dimension of MEG only (gamma1)
S.x7dim = size(V{1,2}{1,2},2);  % Dimension of MEG only (gamma2)

% Classification step
LI = [0.1606 LI(7) LI(8) LI(9) LI(10) LI(11) LI(12)];
LL = [1.0000 LL(7) LL(8) LL(9) LL(10) LL(11) LL(12)];

rng('default') % For reproducibility
[~,acc2] = mkl_ens_n(V,labels,'Nfold',5,'Nrun',1000,'feat_norm',1,'ens',1,...
    'machine','kernel','HyperI',LI,'HyperL',LL,...
    'Hyper_sizeS',S);

save MRIMEG_GrdVar acc2
tbl.VAR_of_GRD = [mean(mean(acc2(:,2:end,:),3))' std(mean(acc2(:,2:end,:),3))'];

% Plot resluts (Classification accuracy and Pos-hoc comparison)

% define titles in same order of the input cell V
titles = {'MRI','MRI,MEG(GRD VAR Delta)','MRI,MEG(GRD VAR Theta)','MRI,MEG(GRD VAR Alpha)',...
    'MRI,MEG(GRD VAR Beta)','MRI,MEG(GRD VAR lGamma)','MRI,MEG(GRD VAR hGamma)'};

pos_titles = {'MRI,MEG(GRD VAR Delta)>MRI','MRI,MEG(GRD VAR Theta)>MRI','MRI,MEG(GRD VAR Alpha)>MRI',...
    'MRI,MEG(GRD VAR Beta)>MRI','MRI,MEG(GRD VAR lGamma)>MRI','MRI,MEG(GRD VAR hGamma)>MRI'};

%  define contrasts
c = [-1 1 0  0 0 0  0;
    -1 0 1  0 0 0  0;
    -1 0 0  1 0 0  0;
    
    -1 0 0  0 1 0  0;
    -1 0 0  0 0 1  0;
    -1 0 0  0 0 0  1];

f6 = plot_results(1,titles,acc2,pos_titles,c);
sgtitle('MEG (GRD VAR) MRI Combinations')
sgt.FontSize = 20;

%% Analysis  (Other MEG - MRI combination), MEG: MAG COV
% Import data and define input cell array
clear
MRI = csvread('derived/ROIdata.csv');
labels = csvread('derived/labels.csv');

V = {{csvread('MAGCOVDELTA.csv')},{csvread('MAGCOVTHETA.csv')},{csvread('MAGCOVALPHA.csv')},...
    {csvread('MAGCOVBETA.csv')},{csvread('MAGCOVGAMMA1.csv')},{csvread('MAGCOVGAMMA2.csv')},...
    {MRI,csvread('MAGCOVDELTA.csv')},{MRI,csvread('MAGCOVTHETA.csv')},{MRI,csvread('MAGCOVALPHA.csv')}...
    ,{MRI,csvread('MAGCOVBETA.csv')},{MRI,csvread('MAGCOVGAMMA1.csv')},{MRI,csvread('MAGCOVGAMMA2.csv')}};

% Tune Hyper-parameters
rng('default') % For reproducibility
[bLcomI,bLcomL,~,~] = mkl_ens_nestedcv(V,labels,'Nfold_o',5,'Nfold_i',10,'Nrun',10);

LI = mean(mean(bLcomI,3)); % Lambda for Intermediate
%  0.5396    0.4232    0.4050    0.3376    0.2656    0.3328    0.4050
%  0.2396    0.1654    0.1996    0.1404    0.1334
LL = mean(mean(bLcomL,3)); % Lambda for Late
%  1.0000    1.0000    1.0000    1.0000    1.0000    1.0000    0.3642
%  0.3790    0.5868    0.6458    0.5372    0.7506

V = {{MRI},{MRI,csvread('MAGCOVDELTA.csv')},{MRI,csvread('MAGCOVTHETA.csv')},{MRI,csvread('MAGCOVALPHA.csv')}...
    ,{MRI,csvread('MAGCOVBETA.csv')},{MRI,csvread('MAGCOVGAMMA1.csv')},{MRI,csvread('MAGCOVGAMMA2.csv')}};

S = struct;
S.x1 = 0.1606; % hyperparameter for MRI only (taken from the first analysis)
S.x2 = LI(1);
S.x3 = LI(2);
S.x4 = LI(3);
S.x5 = LI(4);
S.x6 = LI(5);
S.x7 = LI(6);
S.x1dim = size(MRI,2); % Dimension of MRI only
S.x2dim = size(V{1,2}{1,2},2);  % Dimension of MEG only (delta)
S.x3dim = size(V{1,2}{1,2},2);  % Dimension of MEG only (theta)
S.x4dim = size(V{1,2}{1,2},2);  % Dimension of MEG only (alpha)
S.x5dim = size(V{1,2}{1,2},2);  % Dimension of MEG only (beta)
S.x6dim = size(V{1,2}{1,2},2);  % Dimension of MEG only (gamma1)
S.x7dim = size(V{1,2}{1,2},2);  % Dimension of MEG only (gamma2)

% Classification step
LI = [0.1606 LI(7) LI(8) LI(9) LI(10) LI(11) LI(12)];
LL = [1.0000 LL(7) LL(8) LL(9) LL(10) LL(11) LL(12)];

rng('default') % For reproducibility
[~,acc2] = mkl_ens_n(V,labels,'Nfold',5,'Nrun',1000,'feat_norm',1,'ens',1,...
    'machine','kernel','HyperI',LI,'HyperL',LL,...
    'Hyper_sizeS',S);

save MRIMEG_MagCov acc2
tbl.COV_of_MAG = [mean(mean(acc2(:,2:end,:),3))' std(mean(acc2(:,2:end,:),3))'];

% Plot resluts (Classification accuracy and Pos-hoc comparison)

% define titles in same order of the input cell V
titles = {'MRI','MRI,MEG(MAG COV Delta)','MRI,MEG(MAG COV Theta)','MRI,MEG(MAG COV Alpha)',...
    'MRI,MEG(MAG COV Beta)','MRI,MEG(MAG COV lGamma)','MRI,MEG(MAG COV hGamma)'};

pos_titles = {'MRI,MEG(MAG COV Delta)>MRI','MRI,MEG(MAG COV Theta)>MRI','MRI,MEG(MAG COV Alpha)>MRI',...
    'MRI,MEG(MAG COV Beta)>MRI','MRI,MEG(MAG COV lGamma)>MRI','MRI,MEG(MAG COV hGamma)>MRI'};

%  define contrasts
c = [-1 1 0  0 0 0  0;
    -1 0 1  0 0 0  0;
    -1 0 0  1 0 0  0;
    
    -1 0 0  0 1 0  0;
    -1 0 0  0 0 1  0;
    -1 0 0  0 0 0  1];

f7 = plot_results(1,titles,acc2,pos_titles,c);
sgtitle('MEG (MAG COV) MRI Combinations')
sgt.FontSize = 20;
eval(sprintf('print -f%d -dpng CofMRIMEG.png',f1))
eval(sprintf('print -f%d -dpng CofMRIMEGMEGcon.png',f2))

%% Analysis  (Other MEG - MRI combination), MEG: MAG VAR
% Import data and define input cell array
clear
MRI = csvread('derived/ROIdata.csv');
labels = csvread('derived/labels.csv');

V = {{csvread('MAGVARDELTA.csv')},{csvread('MAGVARTHETA.csv')},{csvread('MAGVARALPHA.csv')},...
    {csvread('MAGCOVBETA.csv')},{csvread('MAGCOVGAMMA1.csv')},{csvread('MAGCOVGAMMA2.csv')},...
    {MRI,csvread('MAGVARDELTA.csv')},{MRI,csvread('MAGVARTHETA.csv')},{MRI,csvread('MAGVARALPHA.csv')}...
    ,{MRI,csvread('MAGVARBETA.csv')},{MRI,csvread('MAGVARGAMMA1.csv')},{MRI,csvread('MAGVARGAMMA2.csv')}};

% Tune Hyper-parameters
rng('default') % For reproducibility
[bLcomI,bLcomL,~,~] = mkl_ens_nestedcv(V,labels,'Nfold_o',5,'Nfold_i',10,'Nrun',10);

LI = mean(mean(bLcomI,3)); % Lambda for Intermediate
% 0.6640    0.6270    0.4702    0.3376    0.2656    0.3328    0.4940
% 0.4844    0.2354    0.2428    0.2132    0.2676
LL = mean(mean(bLcomL,3)); % Lambda for Late
% 1.0000    1.0000    1.0000    1.0000    1.0000    1.0000    0.3826
% 0.3982    0.5862    0.6142    0.5044    0.4478


V = {{MRI},{MRI,csvread('MAGVARDELTA.csv')},{MRI,csvread('MAGVARTHETA.csv')},{MRI,csvread('MAGVARALPHA.csv')}...
    ,{MRI,csvread('MAGVARBETA.csv')},{MRI,csvread('MAGVARGAMMA1.csv')},{MRI,csvread('MAGVARGAMMA2.csv')}};

S = struct;
S.x1 = 0.1606; % hyperparameter for MRI only (taken from the first analysis)
S.x2 = LI(1);
S.x3 = LI(2);
S.x4 = LI(3);
S.x5 = LI(4);
S.x6 = LI(5);
S.x7 = LI(6);
S.x1dim = size(MRI,2); % Dimension of MRI only
S.x2dim = size(V{1,2}{1,2},2);  % Dimension of MEG only (delta)
S.x3dim = size(V{1,2}{1,2},2);  % Dimension of MEG only (theta)
S.x4dim = size(V{1,2}{1,2},2);  % Dimension of MEG only (alpha)
S.x5dim = size(V{1,2}{1,2},2);  % Dimension of MEG only (beta)
S.x6dim = size(V{1,2}{1,2},2);  % Dimension of MEG only (gamma1)
S.x7dim = size(V{1,2}{1,2},2);  % Dimension of MEG only (gamma2)

% Classification step
LI = [0.1606 LI(7) LI(8) LI(9) LI(10) LI(11) LI(12)];
LL = [1.0000 LL(7) LL(8) LL(9) LL(10) LL(11) LL(12)];

% Classification step
rng('default') % For reproducibility
[~,acc2] = mkl_ens_n(V,labels,'Nfold',5,'Nrun',1000,'feat_norm',1,'ens',1,...
    'machine','kernel','HyperI',LI,'HyperL',LL,...
    'Hyper_sizeS',S);

save MRIMEG_MagVar acc2
tbl.VAR_of_MAG = [mean(mean(acc2(:,2:end,:),3))' std(mean(acc2(:,2:end,:),3))'] % table 2 of main paper

% Plot resluts (Classification accuracy and Pos-hoc comparison)

% define titles in same order of the input cell V
titles = {'MRI','MRI,MEG(MAG VAR Delta)','MRI,MEG(MAG VAR Theta)','MRI,MEG(MAG VAR Alpha)',...
    'MRI,MEG(MAG VAR Beta)','MRI,MEG(MAG VAR lGamma)','MRI,MEG(MAG VAR hGamma)'};

pos_titles = {'MRI,MEG(MAG VAR Delta)>MRI','MRI,MEG(MAG VAR Theta)>MRI','MRI,MEG(MAG VAR Alpha)>MRI',...
    'MRI,MEG(MAG VAR Beta)>MRI','MRI,MEG(MAG VAR lGamma)>MRI','MRI,MEG(MAG VAR hGamma)>MRI'};

%  define contrasts
c = [-1 1 0  0 0 0  0;
    -1 0 1  0 0 0  0;
    -1 0 0  1 0 0  0;
    
    -1 0 0  0 1 0  0;
    -1 0 0  0 0 1  0;
    -1 0 0  0 0 0  1];

f8 = plot_results(1,titles,acc2,pos_titles,c);
sgtitle('MEG (MAG VAR) MRI Combinations')
sgt.FontSize = 20;

tbl

%% LOO vs 5-Fold cv for MRI

MRI = csvread('derived/ROIdata.csv');
labels = csvread('derived/labels.csv');

rng(1)
cvp = cvpartition(y,'LeaveOut');
opt = struct("CVPartition",cvp,"ShowPlots",false);
MRImodel = fitcknn(MRI,labels,"NumNeighbors",13,'DistanceWeight', 'squaredinverse','Distance', 'correlation',"CVPartition",cvp);
1-kfoldLoss(MRImodel)

rng(1)
cvp = cvpartition(y,'KFold',5);
opt = struct("CVPartition",cvp,"ShowPlots",false);
MRImodel = fitcknn(MRI,labels,"NumNeighbors",13,'DistanceWeight', 'squaredinverse','Distance', 'correlation',"CVPartition",cvp);
1-kfoldLoss(MRImodel)
