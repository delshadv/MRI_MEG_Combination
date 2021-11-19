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

