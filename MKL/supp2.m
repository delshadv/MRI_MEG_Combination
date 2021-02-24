%% Supplementary results 2
% Combining MEG and MRI using late combination to classify MCI vs Control
% (BioFIND dataset)
%
% This script produces supplementary results after the noise simulations
% paper ""

% Henson R.N 2020, Vaghari D 2020

%% Define necessary paths

bwd = pwd; % needs to be in MRI_MEG_combination directory
addpath(fullfile(bwd,'MKL'));
addpath(fullfile(bwd,'MKL','supplementary'));
addpath(fullfile(bwd,'MKL','derived'));
cd MKL
%% Simple concatenation and Multiple kernel leraning
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

%% Classification accuracy for all MEG features

load (fullfile(bwd,'MKL','derived','labels.csv')) % Classification labels, ie MCI vs CON
tbl = table();
tbl.bands = {'Delta(2-4Hz)';'Theta(4-8Hz)';'Alpha(8-12Hz)';'Beta(12-30Hz)';'Low-Gamma(30-48)';'High-Gamma(52-86)'};

V = {{csvread('GRDCOVDELTA.csv')},{csvread('GRDCOVTHETA.csv')},{csvread('GRDCOVALPHA.csv')}...
    ,{csvread('GRDCOVBETA.csv')},{csvread('GRDCOVGAMMA1.csv')},{csvread('GRDCOVGAMMA2.csv')}};

% Classification step
rng('default') % For reproducibility
[~,acc2] = mkl_ens(V,labels,'Hyper1',0.1,'Hyper2',1,...
    'Nfold',5,'Nrun',1000,'PCA_cut',0,'feat_norm',1,'ens',1);

save MEG_GrdCov acc2
tbl.COV_of_GRD = [mean(mean(acc2,3))' std(mean(acc2,3))'];

clear V
V = {{csvread('GRDVARDELTA.csv')},{csvread('GRDVARTHETA.csv')},{csvread('GRDVARALPHA.csv')}...
    ,{csvread('GRDVARBETA.csv')},{csvread('GRDVARGAMMA1.csv')},{csvread('GRDVARGAMMA2.csv')}};

% Classification step
rng('default') % For reproducibility
[~,acc2] = mkl_ens(V,labels,'Hyper1',0.1,'Hyper2',1,...
    'Nfold',5,'Nrun',1000,'PCA_cut',0,'feat_norm',1,'ens',1);

save MEG_GrdVar acc2
tbl.VAR_of_GRD = [mean(mean(acc2,3))' std(mean(acc2,3))'];

clear V

V = {{csvread('MAGCOVDELTA.csv')},{csvread('MAGCOVTHETA.csv')},{csvread('MAGCOVALPHA.csv')}...
    ,{csvread('MAGCOVBETA.csv')},{csvread('MAGCOVGAMMA1.csv')},{csvread('MAGCOVGAMMA2.csv')}};

% Classification step
rng('default') % For reproducibility
[~,acc2] = mkl_ens(V,labels,'Hyper1',0.1,'Hyper2',1,...
    'Nfold',5,'Nrun',1000,'PCA_cut',0,'feat_norm',1,'ens',1);

save MEG_MagCov acc2
tbl.COV_of_MAG = [mean(mean(acc2,3))' std(mean(acc2,3))'];

clear V

V = {{csvread('MAGVARDELTA.csv')},{csvread('MAGVARTHETA.csv')},{csvread('MAGVARALPHA.csv')}...
    ,{csvread('MAGVARBETA.csv')},{csvread('MAGVARGAMMA1.csv')},{csvread('MAGVARGAMMA2.csv')}};

% Classification step
rng('default') % For reproducibility
[~,acc2] = mkl_ens(V,labels,'Hyper1',0.1,'Hyper2',1,...
    'Nfold',5,'Nrun',1000,'PCA_cut',0,'feat_norm',1,'ens',1);


save MEG_MagVar acc2
tbl.VAR_of_MAG = [mean(mean(acc2,3))' std(mean(acc2,3))']

%% Analysis Supp. (Factorial comparison of MAG vs GRD and VAR vs COV)
% Prepare Input and output for the classifier
clear V
modal = {'MEGMAG','MEGPLANAR'};
load([modal{1} '.mat']) % Needs to run preproc_meg first

V = {};
X1 = variance{5}; % low-gamma
X2 = covariance{5}; % low-gamma

V{1} = {X1};
V{2} = {X2};
V{end+1} = {X1,X2};

load([modal{2} '.mat'])
X1 = variance{5};
X2 = covariance{5};

V{4} = {X1};
V{5} = {X2};
V{end+1} = {X1,X2};

V{7} = {V{1}{1},V{4}{1}}; % VAR MAG,GRD
V{8} = {V{2}{1},V{5}{1}}; % COV MAG,GRD
V{9} = {V{1}{1},V{2}{1},V{4}{1},V{5}{1}}; % VAR,COV MAG,GRD

% Classification step
rng('default') % For reproducibility
[~,acc2] = mkl_ens(V,labels,'Hyper1',0.1,'Hyper2',1,...
    'Nfold',5,'Nrun',1000,'PCA_cut',0,'feat_norm',1,'ens',1);

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
f = plot_results(titles,acc2,pos_titles,c);
sgt = sgtitle('MEG');
sgt.FontSize = 20;

% Save Figures
saveas(f,'Supp Fig5.png')
