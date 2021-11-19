%% Supplementary results 2
% Sections S2-S7 of Supplementaru materials 
% 
% This script produces supplementary results after the noise simulations
% Paper "Late Combination shows that MEG adds to MRI in classifying MCI versus Controls"
% (BioFIND dataset)
%
% Henson R.N 2020, Vaghari D 2020

%% Define necessary paths

bwd = pwd; % needs to be in MRI_MEG_combination directory
addpath(fullfile(bwd,'MKL'));
addpath(fullfile(bwd,'MKL','supplementary'));
addpath(fullfile(bwd,'MKL','derived'));
cd MKL
%% Early and Intermediate Combination (S2)
load(fullfile(bwd,'MKL','CofsMRIMEG_GrdCov')) % needs to run first analysis 1 in main.m 

titles = {'COF','MRI','MEG','MRI,COF','MEG,COF'};
pos_titles = {'MRI>COF','MEG>COF','MEG>MRI',...
            'MRI,COF>COF','MEG,COF>COF'};
%  define contrasts
c = [-1 1 0  0 0 ; 
    -1 0 1   0 0 ; 
    0 -1 1   0 0 ;
    
    -1 0 0   1 0 ;
    -1 0 0   0 1 ];

f2 = plot_results(titles,acc3(:,1:5,:),pos_titles,c); % Supplementary figure 3
sgtitle('Early Combination')
sgt.FontSize = 20;
% Save Figures
saveas(f2,'Supp Fig3.png')

f3 = plot_results(titles,acc1(:,1:5,:),pos_titles,c); % Supplementary figure 4
sgtitle('Intermediate Combination')
sgt.FontSize = 20;

% Save Figures
saveas(f3,'Supp Fig4.png')

%% Classification accuracy for MEG features only (S3)

load (fullfile(bwd,'MKL','derived','labels.csv')) % Classification labels, i.e. MCI vs CON
tbl = table();
tbl.bands = {'Delta(2-4Hz)';'Theta(4-8Hz)';'Alpha(8-12Hz)';'Beta(12-30Hz)';'Low-Gamma(30-48)';'High-Gamma(52-86)'};

clear
labels = csvread('derived/labels.csv');

V = {{csvread('GRDCOVDELTA.csv')},{csvread('GRDCOVTHETA.csv')},{csvread('GRDCOVALPHA.csv')}...
    ,{csvread('GRDCOVBETA.csv')},{csvread('GRDCOVGAMMA1.csv')},{csvread('GRDCOVGAMMA2.csv')}};

LI = [0.1624    0.2892    0.2862    0.3052    0.4490    0.4360]; % taken from previoues analysis
LL = [1 1 1 1 1 1];

S = struct;
S.x1 = LI(1); 
S.x2 = LI(2);
S.x3 = LI(3);
S.x4 = LI(4);
S.x5 = LI(5);
S.x6 = LI(6);
S.x1dim = size(V{1,1}{1,1},2);  % Dimension of MEG only (delta)
S.x2dim = size(V{1,1}{1,1},2);  % Dimension of MEG only (theta)
S.x3dim = size(V{1,1}{1,1},2);  % Dimension of MEG only (alpha)
S.x4dim = size(V{1,1}{1,1},2);  % Dimension of MEG only (beta)
S.x5dim = size(V{1,1}{1,1},2);  % Dimension of MEG only (gamma1)
S.x6dim = size(V{1,1}{1,1},2);  % Dimension of MEG only (gamma2)

% Classification step
rng('default') % For reproducibility
[~,acc2] = mkl_ens_n(V,labels,'Nfold',5,'Nrun',1000,'feat_norm',1,'ens',1,...
    'machine','kernel','HyperI',LI,'HyperL',LL,...
    'Hyper_sizeS',S);

save MEG_GrdCov acc2
tbl.COV_of_GRD = [mean(mean(acc2,3))' std(mean(acc2,3))'];

clear V
V = {{csvread('GRDVARDELTA.csv')},{csvread('GRDVARTHETA.csv')},{csvread('GRDVARALPHA.csv')}...
    ,{csvread('GRDVARBETA.csv')},{csvread('GRDVARGAMMA1.csv')},{csvread('GRDVARGAMMA2.csv')}};

LI = [0.2822    0.3442    0.3168    0.3076    0.2172    0.2968]; % taken from previoues analysis
LL = [1 1 1 1 1 1];

S = struct;
S.x1 = LI(1); 
S.x2 = LI(2);
S.x3 = LI(3);
S.x4 = LI(4);
S.x5 = LI(5);
S.x6 = LI(6);
S.x1dim = size(V{1,1}{1,1},2);  % Dimension of MEG only (delta)
S.x2dim = size(V{1,1}{1,1},2);  % Dimension of MEG only (theta)
S.x3dim = size(V{1,1}{1,1},2);  % Dimension of MEG only (alpha)
S.x4dim = size(V{1,1}{1,1},2);  % Dimension of MEG only (beta)
S.x5dim = size(V{1,1}{1,1},2);  % Dimension of MEG only (gamma1)
S.x6dim = size(V{1,1}{1,1},2);  % Dimension of MEG only (gamma2)

% Classification step
rng('default') % For reproducibility
[~,acc2] = mkl_ens_n(V,labels,'Nfold',5,'Nrun',1000,'feat_norm',1,'ens',1,...
    'machine','kernel','HyperI',LI,'HyperL',LL,...
    'Hyper_sizeS',S);

save MEG_GrdVar acc2
tbl.VAR_of_GRD = [mean(mean(acc2,3))' std(mean(acc2,3))'];

clear V
V = {{csvread('MAGCOVDELTA.csv')},{csvread('MAGCOVTHETA.csv')},{csvread('MAGCOVALPHA.csv')}...
    ,{csvread('MAGCOVBETA.csv')},{csvread('MAGCOVGAMMA1.csv')},{csvread('MAGCOVGAMMA2.csv')}};

LI = [0.5396    0.4232    0.4050    0.3376    0.2656    0.3328]; % taken from previoues analysis
LL = [1 1 1 1 1 1];

S = struct;
S.x1 = LI(1); 
S.x2 = LI(2);
S.x3 = LI(3);
S.x4 = LI(4);
S.x5 = LI(5);
S.x6 = LI(6);
S.x1dim = size(V{1,1}{1,1},2);  % Dimension of MEG only (delta)
S.x2dim = size(V{1,1}{1,1},2);  % Dimension of MEG only (theta)
S.x3dim = size(V{1,1}{1,1},2);  % Dimension of MEG only (alpha)
S.x4dim = size(V{1,1}{1,1},2);  % Dimension of MEG only (beta)
S.x5dim = size(V{1,1}{1,1},2);  % Dimension of MEG only (gamma1)
S.x6dim = size(V{1,1}{1,1},2);  % Dimension of MEG only (gamma2)

% Classification step
rng('default') % For reproducibility
[~,acc2] = mkl_ens_n(V,labels,'Nfold',5,'Nrun',1000,'feat_norm',1,'ens',1,...
    'machine','kernel','HyperI',LI,'HyperL',LL,...
    'Hyper_sizeS',S);

save MEG_MagCov acc2
tbl.COV_of_MAG = [mean(mean(acc2,3))' std(mean(acc2,3))'];

clear V

V = {{csvread('MAGVARDELTA.csv')},{csvread('MAGVARTHETA.csv')},{csvread('MAGVARALPHA.csv')}...
    ,{csvread('MAGVARBETA.csv')},{csvread('MAGVARGAMMA1.csv')},{csvread('MAGVARGAMMA2.csv')}};

LI = [0.6640    0.6270    0.4702    0.3376    0.2656    0.3328]; % taken from previoues analysis
LL = [1 1 1 1 1 1];

S = struct;
S.x1 = LI(1); 
S.x2 = LI(2);
S.x3 = LI(3);
S.x4 = LI(4);
S.x5 = LI(5);
S.x6 = LI(6);
S.x1dim = size(V{1,1}{1,1},2);  % Dimension of MEG only (delta)
S.x2dim = size(V{1,1}{1,1},2);  % Dimension of MEG only (theta)
S.x3dim = size(V{1,1}{1,1},2);  % Dimension of MEG only (alpha)
S.x4dim = size(V{1,1}{1,1},2);  % Dimension of MEG only (beta)
S.x5dim = size(V{1,1}{1,1},2);  % Dimension of MEG only (gamma1)
S.x6dim = size(V{1,1}{1,1},2);  % Dimension of MEG only (gamma2)

% Classification step
rng('default') % For reproducibility
[~,acc2] = mkl_ens_n(V,labels,'Nfold',5,'Nrun',1000,'feat_norm',1,'ens',1,...
    'machine','kernel','HyperI',LI,'HyperL',LL,...
    'Hyper_sizeS',S);

save MEG_MagVar acc2
tbl.VAR_of_MAG = [mean(mean(acc2,3))' std(mean(acc2,3))']

%% Analysis Supp. (Factorial comparison of MAG vs GRD and VAR vs COV)
% Prepare Input and output for the classifier
clear
modal = {'MEGMAG','MEGPLANAR'};
load([modal{1} '.mat'])

V = {};
X1 = variance{6};
X2 = covariance{6};

V{1} = {X1};
V{2} = {X2};
V{end+1} = {X1,X2};

load([modal{2} '.mat'])
X1 = variance{6};
X2 = covariance{6};

V{4} = {X1};
V{5} = {X2};
V{end+1} = {X1,X2};

V{7} = {V{1}{1},V{4}{1}}; % VAR MAG,GRD
V{8} = {V{2}{1},V{5}{1}}; % COV MAG,GRD
V{9} = {V{1}{1},V{2}{1},V{4}{1},V{5}{1}}; % VAR,COV MAG,GRD

load ('derived/labels.csv') % Classification labels, ie MCI vs CON

% Tune Hyper-parameters
rng('default') % For reproducibility
[bLcomI,bLcomL,~,~] = mkl_ens_nestedcv(V,labels,'Nfold_o',5,'Nfold_i',10,'Nrun',10);

LI = mean(mean(bLcomI,3)); % Lambda for Intermediate
% 0.3570    0.3612    0.0212    0.2304    0.4210    0.1190    0.1024    0.0740    0.1334
LL = mean(mean(bLcomL,3)); % Lambda for Late
% 1.0000    1.0000    0.7208    1.0000    1.0000    0.6424    0.5076    0.5196    0.5994

S = struct;
S.x1 = LI(1); % MAG VAR
S.x2 = LI(2); % MAG COV
S.x3 = LI(4); % GRD VAR
S.x4 = LI(5); % GRD COV
S.x1dim = size(V{1,1}{1,1},2);  % Dimension of MEG only (MAG VAR)
S.x2dim = size(V{1,2}{1,1},2);  % Dimension of MEG only (MAG COV)
S.x3dim = size(V{1,4}{1,1},2);  % Dimension of MEG only (GRD VAR)
S.x4dim = size(V{1,5}{1,1},2);  % Dimension of MEG only (GRD COV)

% Classification step
rng('default') % For reproducibility
[~,acc2] = mkl_ens_n(V,labels,'Nfold',5,'Nrun',1000,'feat_norm',1,'ens',1,...
    'machine','kernel','HyperI',LI,'HyperL',LL,...
    'Hyper_sizeS',S);

save MagGrd_VarCov acc2

%titles = {'VAR:GRD','VAR:MAG','VAR:GRD,MAG','COV:GRD','COV:MAG','COV:GRD,MAG','VAR,COV:GRD','VAR,COV:MAG','VAR,COV:GRD,MAG'};
titles = {'VAR:MAG','COV:MAG','VAR,COV:MAG','VAR:GRD','COV:GRD','VAR,COV:GRD','VAR:MAG,GRD','COV:MAG,GRD','VAR,COV:GRD,MAG'};
pos_titles = {'COV>VAR?','COV,VAR>VAR?','GRD>MAG?','GRD,MAG>MAG?'};
%  define contrasts
c = [-1 1 0  -1 1 0  -1 1 0;    % COV > VAR
    -1 0 1  -1 0 1  -1 0 1;   % COV,VAR > VAR
    -1 -1 -1  1 1 1   0 0 0;    % GRD > MAG
    -1 -1 -1  0 0 0  1 1 1;   % GRD,MAG > MAG
    ];
f10 = plot_results(0,titles,acc2,pos_titles,c);
sgt = sgtitle('MEG');
sgt.FontSize = 20;
saveas(f10,'Supp Fig5.png')

%% S4(PCA)
% Import data and define input cell array
MRI = csvread('derived/ROIdata.csv');
MEG = csvread('derived/GRDCOVGAMMA1.csv'); % Since lgamma does best numerically
labels = csvread('derived/labels.csv');
COFS = csvread('derived/COFS.csv');

COFS_cell = {COFS(:,1),COFS(:,2),COFS(:,3),COFS(:,4),COFS(:,4),COFS(:,6),COFS(:,7),COFS(:,8)};

V = {{COFS(:,1)},{MRI},{MEG},{MRI,MEG},[COFS_cell(:);MRI;MEG]',[COFS_cell(:);MRI]'};

% Tune Hyper-parameters
rng('default') % For reproducibility
[bLcomI,bLcomL,~,~] = mkl_ens_nestedcv(V,labels,'Nfold_o',5,'Nfold_i',10,'Nrun',10,'PCA_cut',95);

LI = mean(mean(bLcomI,3)); % Lambda for Intermediate
% 0.4998    0.2000    0.3662    0.1398    0.1436    0.1100
LL = mean(mean(bLcomL,3)); % Lambda for Late
% 1.0000    1.0000    1.0000    0.6228    0.6182    0.4392

S = struct;
S.x1 = LI(1); % hyperparameter for COF only
S.x2 = LI(2); % 0.1606 hyperparameter for MRI only
S.x3 = LI(3); % 0.4000; hyperparameter for MEG only

% Classification step
V = {{MRI},{MEG},{MRI,MEG},[COFS_cell(:);MRI;MEG]',[COFS_cell(:);MRI]'};

LI(1) = [];
LL(1) = [];

rng('default') % For reproducibility
[acc1,acc2] = mkl_ens_n(V,labels,'Nfold',5,'Nrun',1000,'feat_norm',1,'ens',1,...
    'machine','kernel','HyperI',LI,'HyperL',LL,...
    'Hyper_sizeS',S,'PCA_cut',95);

clear V
V = {{[MRI MEG]},{[COFS MRI MEG]}};

% Tune Hyper-parameters
rng('default') % For reproducibility
[bLcomI,~,~,~] = mkl_ens_nestedcv(V,labels,'Nfold_o',5,'Nfold_i',10,'Nrun',10,'ens',0,'PCA_cut',95);

LI = mean(mean(bLcomI,3)); % Lambda for Intermediate
%   0.3754    0.4010

LI = [0.2000 0.3662 LI]; % MRI and MEG (only) taken from the prev analysis 

V = {{MRI},{MEG},{[MRI MEG]},{[COFS MRI MEG]}};

rng('default') % For reproducibility
acc3 = mkl_ens_n(V,labels,'Nfold',5,'Nrun',1000,'feat_norm',1,'ens',0,...
    'machine','kernel','HyperI',LI,'PCA_cut',95);

save PCACofsMRIMEG_GrdCov acc1 acc2 acc3

accuracy = [];

accuracy = cat(2,acc3(:,3,:),acc1(:,3,:));
accuracy = cat(2,accuracy,acc2(:,[3 4 5 1],:));

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
     0 0 1  0 0 -1;
     0 0 0 1 -1 0];

f18 = plot_results(0,titles,accuracy(:,1:end-2,:),pos_titles,c,accuracy);
sgtitle('Early vs Intermediate vs Late Combination (PCA applied)')
sgt.FontSize = 20;
saveas(f18,'Supp PCA.png')

%% S5 (Other classifiers)
% Comparison of MRI vs MEG,MRI using other classifiers
% Import data and define input cell array
MRI = csvread('derived/ROIdata.csv');
MEG = csvread('derived/GRDCOVGAMMA1.csv'); % Since lgamma does best numerically
labels = csvread('derived/labels.csv');

V = {{MRI},{MEG},{MRI,MEG}};

LL = [1 1 0.7416];
% Classification step
rng('default') % For reproducibility
[~,acc_forest] = mkl_ens_n(V,labels,'Nfold',5,'Nrun',1000,'feat_norm',1,'ens',1,...
    'machine','forest','HyperL',LL,'HyperI',LL);
save acc_forest acc_forest

titles = {'MRI','MEG','MRI,MEG'};
pos_titles = {'MRI,MEG>MRI'};
%  define contrasts
c = [-1 0 1];
f11 = plot_results(0,titles,acc_forest,pos_titles,c); % main figure 1
sgt = sgtitle('Random Forest');
sgt.FontSize = 20;
% Save Figures
saveas(f11,'Supp Fig_forest.png')
% 0.9290 (3.9696)
% 68.2438   65.7121   72.2133

% Classification step
rng('default') % For reproducibility
[~,acc_knn] = mkl_ens_n(V,labels,'Nfold',5,'Nrun',1000,'feat_norm',1,'ens',1,...
    'machine','knn','HyperL',LL,'HyperI',LL);
save acc_knn acc_knn

f12 = plot_results(0,titles,acc_knn,pos_titles,c); % main figure 1
sgt = sgtitle('KNN');
sgt.FontSize = 20;
% Save Figures
saveas(f12,'Supp Fig_KNN.png')

rng('default') % For reproducibility
[~,acc_net] = mkl_ens_n(V,labels,'Nfold',5,'Nrun',1000,'feat_norm',1,'ens',1,...
    'machine','deepNet','HyperL',LL,'HyperI',LL);
save acc_net acc_net

f13 = plot_results(0,titles,acc_net,pos_titles,c); % main figure 1
sgt = sgtitle('Neural Net');
sgt.FontSize = 20;
% Save Figures
saveas(f13,'Supp FigDNN.png')

%% S6 (MRI Voxels)
%% Analysis Supp (Comparison of MRI vs MEG when MRI contains all voxels)
% Import data and define input cell array
MRI = csvread('derived/MRIVOXELS.csv');
MEG = csvread('derived/GRDCOVGAMMA1.csv'); % Since lgamma does best numerically
labels = csvread('derived/labels.csv');

V = {{MRI},{MRI,MEG}};

% Tune Hyper-parameters
rng('default') % For reproducibility
[bLcomI,bLcomL,~,~] = mkl_ens_nestedcv(V,labels,'Nfold_o',5,'Nfold_i',10,'Nrun',10);

LI = mean(mean(bLcomI,3)); % Lambda for Intermediate
% 0.2192   0.1350
LL = mean(mean(bLcomL,3)); % Lambda for Late
% 1.0000   0.5530

S = struct;
S.x1 = LI(1); % hyperparameter for MRI only
S.x1dim = size(MRI,2);  % Dimension of MRI only
S.x2 = 0.4;   % hyperparameter for MEG only (taken from the first analysis)
S.x2dim = size(MEG,2); % Dimension of MEG only

V = {{MRI},{MEG},{MRI,MEG}};

LI = [LI(1) 0.4 LI(2)];
LL = [LL(1) 1.0 LL(2)];

% Classification step
rng('default') % For reproducibility
[acc1,acc2] = mkl_ens_n(V,labels,'Nfold',5,'Nrun',10,'feat_norm',1,'ens',1,...
    'machine','kernel','HyperI',LI,'HyperL',LL,...
    'Hyper_sizeS',S);

save acc_voxel acc1 acc2

titles = {'MRI','MEG','MRI,MEG'};
pos_titles = {'MRI,MEG>MRI'};

%  define contrasts
c = [-1 0 1];
f14 = plot_results(0,titles,acc1,pos_titles,c); % main figure 1
sgt = sgtitle('MEG,MRI (Voxels)');
sgt.FontSize = 20;
% Save Figures
saveas(f14,'Supp Fig_voxels.png')

%% S7 (converter vs non-converter MCI)
%% Converters vs Non-Converters (GRD-COV)

clear
MRI = csvread('derived/ROIdata.csv');
conv_indx = csvread('derived/conv_indx.csv');

MRI(isnan(conv_indx),:) = [];

MEG1 = csvread('derived/GRDCOVDELTA.csv'); MEG1(isnan(conv_indx),:) = [];
MEG2 = csvread('derived/GRDCOVTHETA.csv'); MEG2(isnan(conv_indx),:) = [];
MEG3 = csvread('derived/GRDCOVALPHA.csv'); MEG3(isnan(conv_indx),:) = [];
MEG4 = csvread('derived/GRDCOVBETA.csv'); MEG4(isnan(conv_indx),:) = [];
MEG5 = csvread('derived/GRDCOVGAMMA1.csv'); MEG5(isnan(conv_indx),:) = [];
MEG6 = csvread('derived/GRDCOVGAMMA2.csv'); MEG6(isnan(conv_indx),:) = [];

conv = conv_indx;
conv(isnan(conv_indx),:) = [];
labels = conv_indx;
labels(isnan(conv_indx),:) = [];
labels = conv+1; % 1 means stable

V = {{MRI},{MEG1},{MEG2},{MEG3},{MEG4},{MEG5},{MEG6},...
    {MRI,MEG1},{MRI,MEG2},{MRI,MEG3}...
    ,{MRI,MEG4},{MRI,MEG5},{MRI,MEG6}};

% Tune Hyper-parameters
rng('default') % For reproducibility
[bLcomI,bLcomL,~,~] = mkl_ens_nestedcv(V,labels,'Nfold_o',5,'Nfold_i',10,'Nrun',10);

LI = mean(mean(bLcomI,3)); % Lambda for Intermediate
% 0.3794    0.4524    0.3256    0.6134    0.6040    0.6200    0.5962
% 0.1486    0.0734    0.1230    0.0212    0.0684    0.0156
LL = mean(mean(bLcomL,3)); % Lambda for Late
%  1.0000    1.0000    1.0000    1.0000    1.0000    1.0000    1.0000
%  0.5246    0.5804    0.6726    0.3464    0.4884    0.2982 

V = {{MRI},{MRI,MEG1},{MRI,MEG2},{MRI,MEG3}...
    ,{MRI,MEG4},{MRI,MEG5},{MRI,MEG6}};

S = struct;
S.x1 = LI(1); % hyperparameter for MRI only 
S.x2 = LI(2);
S.x3 = LI(3);
S.x4 = LI(4);
S.x5 = LI(5);
S.x6 = LI(6);
S.x7 = LI(7);
S.x1dim = size(MRI,2); % Dimension of MRI only
S.x2dim = size(MEG1,2);  % Dimension of MEG only (delta)
S.x3dim = size(MEG1,2);  % Dimension of MEG only (theta)
S.x4dim = size(MEG1,2);  % Dimension of MEG only (alpha)
S.x5dim = size(MEG1,2);  % Dimension of MEG only (beta)
S.x6dim = size(MEG1,2);  % Dimension of MEG only (gamma1)
S.x7dim = size(MEG1,2);  % Dimension of MEG only (gamma2)

% Classification step
LI (1:7) = []; LI = [0.3794 LI];
LL (1:7) = []; LL = [1.0000 LL];

rng('default') % For reproducibility
[~,acc2] = mkl_ens_n(V,labels,'Nfold',5,'Nrun',1000,'feat_norm',1,'ens',1,...
    'machine','kernel','HyperI',LI,'HyperL',LL,...
    'Hyper_sizeS',S);

save ConvMRIMEG_GrdCov acc2
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

f14 = plot_results(1,titles,acc2,pos_titles,c);
sgtitle('MEG (GRD COV) MRI Combinations, Converter')
sgt.FontSize = 20;

% 0.8230    0.7450    0.6690    0.5800    0.5180    0.5760
% 3.0493    2.6783    1.5189    0.7741    0.3552    0.8665
% 68.7119   71.7612   71.3902   70.2308   69.4860   69.0671   69.5784 Mean
% 4.9984    4.9637    4.8482    5.0413    5.1183    4.8461    4.8522  Std
%% Converters vs Non-Converters (MAG-COV)

clear
MRI = csvread('derived/ROIdata.csv');
conv_indx = csvread('derived/conv_indx.csv');

MRI(isnan(conv_indx),:) = [];

MEG1 = csvread('derived/MAGCOVDELTA.csv'); MEG1(isnan(conv_indx),:) = [];
MEG2 = csvread('derived/MAGCOVTHETA.csv'); MEG2(isnan(conv_indx),:) = [];
MEG3 = csvread('derived/MAGCOVALPHA.csv'); MEG3(isnan(conv_indx),:) = [];
MEG4 = csvread('derived/MAGCOVBETA.csv'); MEG4(isnan(conv_indx),:) = [];
MEG5 = csvread('derived/MAGCOVGAMMA1.csv'); MEG5(isnan(conv_indx),:) = [];
MEG6 = csvread('derived/MAGCOVGAMMA2.csv'); MEG6(isnan(conv_indx),:) = [];

conv = conv_indx;
conv(isnan(conv_indx),:) = [];
labels = conv_indx;
labels(isnan(conv_indx),:) = [];
labels = conv+1; % 1 means stable

V = {{MRI},{MEG1},{MEG2},{MEG3},{MEG4},{MEG5},{MEG6},...
    {MRI,MEG1},{MRI,MEG2},{MRI,MEG3}...
    ,{MRI,MEG4},{MRI,MEG5},{MRI,MEG6}};

% Tune Hyper-parameters
rng('default') % For reproducibility
[bLcomI,bLcomL,~,~] = mkl_ens_nestedcv(V,labels,'Nfold_o',5,'Nfold_i',10,'Nrun',10);

LI = mean(mean(bLcomI,3)); % Lambda for Intermediate
% 0.3334    0.4452    0.2390    0.5644    0.5380    0.4642    0.4840
% 0.1158    0.1452    0.0726    0.0946    0.0736    0.1348
LL = mean(mean(bLcomL,3)); % Lambda for Late
% 1.0000    1.0000    1.0000    1.0000    1.0000    1.0000    1.0000
% 0.5190    0.3682    0.6330    0.4026    0.6370    0.5318

V = {{MRI},{MRI,MEG1},{MRI,MEG2},{MRI,MEG3}...
    ,{MRI,MEG4},{MRI,MEG5},{MRI,MEG6}};

S = struct;
S.x1 = LI(1); % hyperparameter for MRI only 
S.x2 = LI(2);
S.x3 = LI(3);
S.x4 = LI(4);
S.x5 = LI(5);
S.x6 = LI(6);
S.x7 = LI(7);
S.x1dim = size(MRI,2); % Dimension of MRI only
S.x2dim = size(MEG1,2);  % Dimension of MEG only (delta)
S.x3dim = size(MEG1,2);  % Dimension of MEG only (theta)
S.x4dim = size(MEG1,2);  % Dimension of MEG only (alpha)
S.x5dim = size(MEG1,2);  % Dimension of MEG only (beta)
S.x6dim = size(MEG1,2);  % Dimension of MEG only (gamma1)
S.x7dim = size(MEG1,2);  % Dimension of MEG only (gamma2)

% Classification step
LI (1:7) = []; LI = [0.3334 LI];
LL (1:7) = []; LL = [1.0000 LL];

rng('default') % For reproducibility
[~,acc2] = mkl_ens_n(V,labels,'Nfold',5,'Nrun',1000,'feat_norm',1,'ens',1,...
    'machine','kernel','HyperI',LI,'HyperL',LL,...
    'Hyper_sizeS',S);

save ConvMRIMEG_MagCov acc2
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

f15 = plot_results(1,titles,acc2,pos_titles,c);
sgtitle('MEG (MAG COV) MRI Combinations, Converter')
sgt.FontSize = 20;

% 0.8220    0.7490    0.6880    0.6420    0.4630    0.6710
% 2.1770    2.4781    1.1775    0.7142   -0.3252    1.4219
% 68.7713   70.9484   71.2495   69.9488   69.4855   68.4461   70.1933    Mean
% 4.7536    4.7805    4.6398    4.8695    4.9297    4.7451    4.5804     Std
%% Converters vs Non-Converters (GRD-VAR)

clear
MRI = csvread('derived/ROIdata.csv');
conv_indx = csvread('derived/conv_indx.csv');

MRI(isnan(conv_indx),:) = [];

MEG1 = csvread('derived/GRDVARDELTA.csv'); MEG1(isnan(conv_indx),:) = [];
MEG2 = csvread('derived/GRDVARTHETA.csv'); MEG2(isnan(conv_indx),:) = [];
MEG3 = csvread('derived/GRDVARALPHA.csv'); MEG3(isnan(conv_indx),:) = [];
MEG4 = csvread('derived/GRDVARBETA.csv'); MEG4(isnan(conv_indx),:) = [];
MEG5 = csvread('derived/GRDVARGAMMA1.csv'); MEG5(isnan(conv_indx),:) = [];
MEG6 = csvread('derived/GRDVARGAMMA2.csv'); MEG6(isnan(conv_indx),:) = [];

conv = conv_indx;
conv(isnan(conv_indx),:) = [];
labels = conv_indx;
labels(isnan(conv_indx),:) = [];
labels = conv+1; % 1 means stable

V = {{MRI},{MEG1},{MEG2},{MEG3},{MEG4},{MEG5},{MEG6},...
    {MRI,MEG1},{MRI,MEG2},{MRI,MEG3}...
    ,{MRI,MEG4},{MRI,MEG5},{MRI,MEG6}};

% Tune Hyper-parameters
rng('default') % For reproducibility
[bLcomI,bLcomL,~,~] = mkl_ens_nestedcv(V,labels,'Nfold_o',5,'Nfold_i',10,'Nrun',10);

LI = mean(mean(bLcomI,3)); % Lambda for Intermediate
% 0.3794    0.2686    0.4254    0.6092    0.5210    0.5780    0.5020
% 0.1388    0.1624    0.1318    0.0732    0.0466    0.0776
LL = mean(mean(bLcomL,3)); % Lambda for Late
% 1.0000    1.0000    1.0000    1.0000    1.0000    1.0000    1.0000
% 0.3788    0.3382    0.5284    0.3428    0.3460    0.3032

V = {{MRI},{MRI,MEG1},{MRI,MEG2},{MRI,MEG3}...
    ,{MRI,MEG4},{MRI,MEG5},{MRI,MEG6}};

S = struct;
S.x1 = LI(1); % hyperparameter for MRI only 
S.x2 = LI(2);
S.x3 = LI(3);
S.x4 = LI(4);
S.x5 = LI(5);
S.x6 = LI(6);
S.x7 = LI(7);
S.x1dim = size(MRI,2); % Dimension of MRI only
S.x2dim = size(MEG1,2);  % Dimension of MEG only (delta)
S.x3dim = size(MEG1,2);  % Dimension of MEG only (theta)
S.x4dim = size(MEG1,2);  % Dimension of MEG only (alpha)
S.x5dim = size(MEG1,2);  % Dimension of MEG only (beta)
S.x6dim = size(MEG1,2);  % Dimension of MEG only (gamma1)
S.x7dim = size(MEG1,2);  % Dimension of MEG only (gamma2)

% Classification step
LI (1:7) = []; LI = [0.3794 LI];
LL (1:7) = []; LL = [1.0000 LL];

rng('default') % For reproducibility
[~,acc2] = mkl_ens_n(V,labels,'Nfold',5,'Nrun',1000,'feat_norm',1,'ens',1,...
    'machine','kernel','HyperI',LI,'HyperL',LL,...
    'Hyper_sizeS',S);

save ConvMRIMEG_GrdVar acc2
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

f16 = plot_results(1,titles,acc2,pos_titles,c);
sgtitle('MEG (GRD VAR) MRI Combinations, Converter')
sgt.FontSize = 20;


% 0.6250    0.6200    0.6540    0.5940    0.6510    0.6500
% 0.9097    0.8089    0.9547    0.6999    1.1464    1.5510
% 68.7499   69.6596   69.5588   69.7046   69.4498   69.8963   70.3009    Mean
% 5.0095    4.9670    4.9403    5.0873    5.0294    4.9038    5.0538     Std
%% Converters vs Non-Converters (MAG-VAR)

clear
MRI = csvread('derived/ROIdata.csv');
conv_indx = csvread('derived/conv_indx.csv');

MRI(isnan(conv_indx),:) = [];

MEG1 = csvread('derived/MAGVARDELTA.csv'); MEG1(isnan(conv_indx),:) = [];
MEG2 = csvread('derived/MAGVARTHETA.csv'); MEG2(isnan(conv_indx),:) = [];
MEG3 = csvread('derived/MAGVARALPHA.csv'); MEG3(isnan(conv_indx),:) = [];
MEG4 = csvread('derived/MAGVARBETA.csv'); MEG4(isnan(conv_indx),:) = [];
MEG5 = csvread('derived/MAGVARGAMMA1.csv'); MEG5(isnan(conv_indx),:) = [];
MEG6 = csvread('derived/MAGVARGAMMA2.csv'); MEG6(isnan(conv_indx),:) = [];

conv = conv_indx;
conv(isnan(conv_indx),:) = [];
labels = conv_indx;
labels(isnan(conv_indx),:) = [];
labels = conv+1; % 1 means stable

V = {{MRI},{MEG1},{MEG2},{MEG3},{MEG4},{MEG5},{MEG6},...
    {MRI,MEG1},{MRI,MEG2},{MRI,MEG3}...
    ,{MRI,MEG4},{MRI,MEG5},{MRI,MEG6}};

% Tune Hyper-parameters
rng('default') % For reproducibility
[bLcomI,bLcomL,~,~] = mkl_ens_nestedcv(V,labels,'Nfold_o',5,'Nfold_i',10,'Nrun',10);

LI = mean(mean(bLcomI,3)); % Lambda for Intermediate
% 0.3794    0.2792    0.3644    0.4786    0.5194    0.5870    0.4632
% 0.2128    0.1634    0.0790    0.1166    0.0632    0.0458 
LL = mean(mean(bLcomL,3)); % Lambda for Late
% 1.0000    1.0000    1.0000    1.0000    1.0000    1.0000    1.0000
% 0.3296    0.2628    0.4438    0.3642    0.2968    0.3002

V = {{MRI},{MRI,MEG1},{MRI,MEG2},{MRI,MEG3}...
    ,{MRI,MEG4},{MRI,MEG5},{MRI,MEG6}};

S = struct;
S.x1 = LI(1); % hyperparameter for MRI only 
S.x2 = LI(2);
S.x3 = LI(3);
S.x4 = LI(4);
S.x5 = LI(5);
S.x6 = LI(6);
S.x7 = LI(7);
S.x1dim = size(MRI,2); % Dimension of MRI only
S.x2dim = size(MEG1,2);  % Dimension of MEG only (delta)
S.x3dim = size(MEG1,2);  % Dimension of MEG only (theta)
S.x4dim = size(MEG1,2);  % Dimension of MEG only (alpha)
S.x5dim = size(MEG1,2);  % Dimension of MEG only (beta)
S.x6dim = size(MEG1,2);  % Dimension of MEG only (gamma1)
S.x7dim = size(MEG1,2);  % Dimension of MEG only (gamma2)

% Classification step
LI (1:7) = []; LI = [0.3794 LI];
LL (1:7) = []; LL = [1.0000 LL];

rng('default') % For reproducibility
[~,acc2] = mkl_ens_n(V,labels,'Nfold',5,'Nrun',1000,'feat_norm',1,'ens',1,...
    'machine','kernel','HyperI',LI,'HyperL',LL,...
    'Hyper_sizeS',S);

save ConvMRIMEG_MagVar acc2
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

f17 = plot_results(1,titles,acc2,pos_titles,c);
sgtitle('MEG (MAG VAR) MRI Combinations, Converter')
sgt.FontSize = 20;
