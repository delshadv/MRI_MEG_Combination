% Main Script:
% This is the "main" script to reproduce results of "" published in ""
% This script carries out permutation test for MCI vs Control classification
% problem using late combination of MEG - MRI.

% You might need to use Parallel Computing Toolbox to be able to run MKL functions
% in parallel for many repetitions to get more accurate results.
% 
% by Rik Henson and Delshad Vaghari, 2020

%% Define Path to needed functions

% Assumed you are currently in the directory including BioFIND data,
% OSL and MKL directories as described in readme.md
bwd = pwd; % needs to be in Github directory
addpath(fullfile(bwd,'MKL')); 
addpath(fullfile(bwd,'MKL','supplementary'));
addpath(fullfile(bwd,'MKL','derived'));

cd MKL % to save variables in MKL dir

%% Analysis 1 (Comparing Early vs Intermediate vs Late Combs)

% Prepare input and output of the classifier
%% Analysis 1 (Comparison of MRI vs MEG vs COFs through all possible combination methods)
% Import data and define input cell array
MRI = csvread(fullfile('derived','ROIdata.csv'));
MEG = csvread(fullfile('derived','GRDCOVGAMMA1.csv')); % Since lgamma does best numerically
COFS = csvread(fullfile('derived','COFS.csv'));
labels = csvread(fullfile('derived','labels.csv'));

COFS_cell = {COFS(:,1),COFS(:,2),COFS(:,3),COFS(:,4),COFS(:,4),COFS(:,6),COFS(:,7),COFS(:,8)};
V = {COFS_cell,{MRI},{MEG},[COFS_cell(:);MRI]',[COFS_cell(:);{MEG}]',{MRI,MEG},...
    [COFS_cell(:);MRI;MEG]'};

% Classification step
rng('default') % For reproducibility
[acc1,acc2] = mkl_ens(V,labels,'Hyper1',0.1,'Hyper2',1,...
    'CVratio',[0.8 0.2],'Nrun',1000,'PCA_cut',0,'feat_norm',1,'ens',1);

clear V
rng('default') % For reproducibility
V = {{COFS},{MRI},{MEG},{[COFS MRI]},{[COFS MEG]},{[MRI MEG]},{[COFS MRI MEG]}};
acc3 = mkl_ens(V,labels,'Hyper1',0.1,'Hyper2',1,...
    'CVratio',[0.8 0.2],'Nrun',1000,'PCA_cut',0,'feat_norm',1,'ens',0);

save CofsMRIMEG_GrdCov acc1 acc2 acc3

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

f1 = plot_results(titles,acc2(:,1:5,:),pos_titles,c); % main figure 1
sgtitle('Late Combination')
sgt.FontSize = 20;
saveas(f1,'Fig1.png')
%f2 = plot_results(titles,acc3(:,1:5,:),pos_titles,c); % Supplementary figure 3
%sgtitle('Early Combination')
%sgt.FontSize = 20;
%f3 = plot_results(titles,acc1(:,1:5,:),pos_titles,c); % Supplementary figure 4
%sgtitle('Intermediate Combination')
%sgt.FontSize = 20;

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
 
f4 = plot_results(titles,accuracy(:,1:end-2,:),pos_titles,c,accuracy); % Main Figure 2
sgtitle('Early vs Intermediate vs Late Combination')
sgt.FontSize = 20;
saveas(f4,'Fig2.png')

