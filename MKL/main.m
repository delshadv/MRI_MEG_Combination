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
MRI = csvread('derived/ROIdata.csv');
MEG = csvread('derived/GRDCOVGAMMA1.csv'); % Since lgamma does best numerically
COFS = csvread('derived/COFS.csv');
labels = csvread('derived/labels.csv');

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
f2 = plot_results(titles,acc3(:,1:5,:),pos_titles,c); % Supplementary figure 3
sgtitle('Early Combination')
sgt.FontSize = 20;
f3 = plot_results(titles,acc1(:,1:5,:),pos_titles,c); % Supplementary figure 4
sgtitle('Intermediate Combination')
sgt.FontSize = 20;

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

%% Main Table 2
% Analysis  (Other MEG - MRI combination), MEG: GRD COV
% Import data and define input cell array
clear
MRI = csvread('derived/ROIdata.csv');
labels = csvread('derived/labels.csv');

V = {{MRI},{MRI,csvread('GRDCOVDELTA.csv')},{MRI,csvread('GRDCOVTHETA.csv')},{MRI,csvread('GRDCOVALPHA.csv')}...
    ,{MRI,csvread('GRDCOVBETA.csv')},{MRI,csvread('GRDCOVGAMMA1.csv')},{MRI,csvread('GRDCOVGAMMA2.csv')}};

% Classification step
rng('default') % For reproducibility
[~,acc2] = mkl_ens(V,labels,'Hyper1',0.1,'Hyper2',1,...
    'CVratio',[0.8 0.2],'Nrun',1000,'PCA_cut',0,'feat_norm',1,'ens',1);

save MRIMEG_GrdCov acc2
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

f5 = plot_results(titles,acc2,pos_titles,c);
sgtitle('MEG (GRD COV) MRI Combinations')
sgt.FontSize = 20;

% Analysis  (Other MEG - MRI combination), MEG: GRD VAR
% Import data and define input cell array
clear
MRI = csvread('derived/ROIdata.csv');
labels = csvread('derived/labels.csv');

V = {{MRI},{MRI,csvread('GRDVARDELTA.csv')},{MRI,csvread('GRDVARTHETA.csv')},{MRI,csvread('GRDVARALPHA.csv')}...
    ,{MRI,csvread('GRDVARBETA.csv')},{MRI,csvread('GRDVARGAMMA1.csv')},{MRI,csvread('GRDVARGAMMA2.csv')}};

% Classification step
rng('default') % For reproducibility
[~,acc2] = mkl_ens(V,labels,'Hyper1',0.1,'Hyper2',1,...
    'CVratio',[0.8 0.2],'Nrun',1000,'PCA_cut',0,'feat_norm',1,'ens',1);

save MRIMEG_GrdVar acc2
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

f6 = plot_results(titles,acc2,pos_titles,c);
sgtitle('MEG (GRD VAR) MRI Combinations')
sgt.FontSize = 20;

% Analysis  (Other MEG - MRI combination), MEG: MAG COV
% Import data and define input cell array
clear
MRI = csvread('derived/ROIdata.csv');
labels = csvread('derived/labels.csv');

V = {{MRI},{MRI,csvread('MAGCOVDELTA.csv')},{MRI,csvread('MAGCOVTHETA.csv')},{MRI,csvread('MAGCOVALPHA.csv')}...
    ,{MRI,csvread('MAGCOVBETA.csv')},{MRI,csvread('MAGCOVGAMMA1.csv')},{MRI,csvread('MAGCOVGAMMA2.csv')}};

% Classification step
rng('default') % For reproducibility
[~,acc2] = mkl_ens(V,labels,'Hyper1',0.1,'Hyper2',1,...
    'CVratio',[0.8 0.2],'Nrun',1000,'PCA_cut',0,'feat_norm',1,'ens',1);

save MRIMEG_MagCov acc2
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

f7 = plot_results(titles,acc2,pos_titles,c);
sgtitle('MEG (MAG COV) MRI Combinations')
sgt.FontSize = 20;

% Analysis  (Other MEG - MRI combination), MEG: MAG VAR
% Import data and define input cell array
clear
MRI = csvread('derived/ROIdata.csv');
labels = csvread('derived/labels.csv');

V = {{MRI},{MRI,csvread('MAGVARDELTA.csv')},{MRI,csvread('MAGVARTHETA.csv')},{MRI,csvread('MAGVARALPHA.csv')}...
    ,{MRI,csvread('MAGVARBETA.csv')},{MRI,csvread('MAGVARGAMMA1.csv')},{MRI,csvread('MAGVARGAMMA2.csv')}};

% Classification step
rng('default') % For reproducibility
[~,acc2] = mkl_ens(V,labels,'Hyper1',0.1,'Hyper2',1,...
    'CVratio',[0.8 0.2],'Nrun',1000,'PCA_cut',0,'feat_norm',1,'ens',1);

save MRIMEG_MagVar acc2
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

f8 = plot_results(titles,acc2,pos_titles,c);
sgtitle('MEG (MAG VAR) MRI Combinations')
sgt.FontSize = 20;
eval(sprintf('print -f%d -dpng CofMRIMEG.png',f1))
eval(sprintf('print -f%d -dpng CofMRIMEGMEGcon.png',f2))










%% Supplementary Table 1 
% Supplementary Analyses (Other MEG - MRI combination), MEG: GRD COV
% Import data and define input cell array
clear
labels = csvread('derived/labels.csv');

V = {{csvread('GRDCOVDELTA.csv')},{csvread('GRDCOVTHETA.csv')},{csvread('GRDCOVALPHA.csv')}...
    ,{csvread('GRDCOVBETA.csv')},{csvread('GRDCOVGAMMA1.csv')},{csvread('GRDCOVGAMMA2.csv')}};

% Classification step
rng('default') % For reproducibility
[~,acc2] = mkl_ens(V,labels,'Hyper1',0.1,'Hyper2',1,...
    'CVratio',[0.8 0.2],'Nrun',1000,'PCA_cut',0,'feat_norm',1,'ens',1);

save MEG_GrdCov acc2
mean(mean(acc2,3))
std(mean(acc2,3))

clear V

V = {{csvread('GRDVARDELTA.csv')},{csvread('GRDVARTHETA.csv')},{csvread('GRDVARALPHA.csv')}...
    ,{csvread('GRDVARBETA.csv')},{csvread('GRDVARGAMMA1.csv')},{csvread('GRDVARGAMMA2.csv')}};

% Classification step
rng('default') % For reproducibility
[~,acc2] = mkl_ens(V,labels,'Hyper1',0.1,'Hyper2',1,...
    'CVratio',[0.8 0.2],'Nrun',1000,'PCA_cut',0,'feat_norm',1,'ens',1);

save MEG_GrdVar acc2
mean(mean(acc2,3))
std(mean(acc2,3))

clear V

V = {{csvread('MAGCOVDELTA.csv')},{csvread('MAGCOVTHETA.csv')},{csvread('MAGCOVALPHA.csv')}...
    ,{csvread('MAGCOVBETA.csv')},{csvread('MAGCOVGAMMA1.csv')},{csvread('MAGCOVGAMMA2.csv')}};

% Classification step
rng('default') % For reproducibility
[~,acc2] = mkl_ens(V,labels,'Hyper1',0.1,'Hyper2',1,...
    'CVratio',[0.8 0.2],'Nrun',1000,'PCA_cut',0,'feat_norm',1,'ens',1);

save MEG_MagCov acc2
mean(mean(acc2,3))
std(mean(acc2,3))

clear V

V = {{csvread('MAGVARDELTA.csv')},{csvread('MAGVARTHETA.csv')},{csvread('MAGVARALPHA.csv')}...
    ,{csvread('MAGVARBETA.csv')},{csvread('MAGVARGAMMA1.csv')},{csvread('MAGVARGAMMA2.csv')}};

% Classification step
rng('default') % For reproducibility
[~,acc2] = mkl_ens(V,labels,'Hyper1',0.1,'Hyper2',1,...
    'CVratio',[0.8 0.2],'Nrun',1000,'PCA_cut',0,'feat_norm',1,'ens',1);


save MEG_MagVar acc2
mean(mean(acc2,3))
std(mean(acc2,3))

%% Supplementary Analyses (Factorial comparison of MAG vs GRD and VAR vs COV) 
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

% Classification step
rng('default') % For reproducibility
[~,acc2] = mkl_ens(V,labels,'Hyper1',0.1,'Hyper2',1,...
    'CVratio',[0.8 0.2],'Nrun',1000,'PCA_cut',0,'feat_norm',1,'ens',1);

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
f1 = plot_results(titles,acc2,pos_titles,c); % Supplementary Figure 5
sgt = sgtitle('MEG');
sgt.FontSize = 20;
