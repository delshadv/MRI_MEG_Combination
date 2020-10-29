% Main Script:
% This is the "main" script to reproduce results of "" published in ""
% This script carries out permutation test for MCI vs Control classification
% problem using combination of MEG - MRI  with Multi-Kernel Learning.

% You might need to use Parallel Computing Toolbox to be able to run MKL functions
% in parallel for many repetitions to get more accurate results.
% 
% by Rik Henson and Delshad Vaghari, 2020

% TO DO:    re-run on single broadband power

%% Define Path to needed functions

% Assumed you are currently in the directory including BioFIND data,
% OSL and MKL directories as described in readme.md
bwd = pwd; % needs to be in Github directory
addpath(fullfile(bwd,'MKL')); 
addpath(fullfile(bwd,'MKL','supplementary'));
cd MKL % to save variables in MKL dir

%% Analysis 1 (Factorial comparison of MAG vs GRD and VAR vs COV)

% Prepare input and output of the classifier
clear
modal = {'MEGMAG','MEGPLANAR'};
V = {};
for m = 1:length(modal)
    load([modal{m} '.mat'])
    X1 = [];
    for f = 2:length(variance) % changed to 2 to bypass concatenation element
        X1 = [X1  variance{f}/std(variance{f}(:))];
    end
    V{end+1} = {X1};
    X2 = [];
    for f = 2:length(covariance) % changed to 2 to bypass concatenation element
        X2 = [X2  covariance{f}/std(covariance{f}(:))];
    end
    V{end+1} = {X2};
    V{end+1} = {X1,X2}; 
end
V{7} = {V{1}{1},V{4}{1}}; % VAR MAG,GRD
V{8} = {V{2}{1},V{5}{1}}; % COV MAG,GRD
V{9} = {V{1}{1},V{2}{1},V{4}{1},V{5}{1}}; % VAR,COV MAG,GRD

load ('y.mat') % Classification labels, ie MCI vs CON

% Classification step
rng('default') % For reproducibility
acc = mkl_class(V,y,'machine','easy',...
    'hyper',1,'CVratio',[0.8 0.2],...
    'Nrun',1000,'PCA_cut',0,'norm',1);

save MagGrd_VarCov acc

%titles = {'VAR:GRD','VAR:MAG','VAR:GRD,MAG','COV:GRD','COV:MAG','COV:GRD,MAG','VAR,COV:GRD','VAR,COV:MAG','VAR,COV:GRD,MAG'};
titles = {'VAR:MAG','COV:MAG','VAR,COV:MAG','VAR:GRD','COV:GRD','VAR,COV:GRD','VAR:MAG,GRD','COV:MAG,GRD','VAR,COV:MAG,GRD'};
pos_titles = {'COV>VAR?','COV,VAR>COV?','GRD>MAG?','GRD,MAG>GRD?'};
%  define contrasts
c = [-1 1 0  -1 1 0  -1 1 0;    % COV > VAR
    0 -1 1    0 -1 1  0 -1 1;   % COV,VAR > COV
    -1 -1 -1  1 1 1   0 0 0;    % GRD > MAG
    0 0 0    -1 -1 -1  1 1 1;   % GRD,MAG > GRD
    ];
    
[f1,f2] = plot_results(titles,acc,pos_titles,c);
sgt = sgtitle('MEG','Color',[0.9290 0.6940 0.1250],'FixedWidth'); 
sgt.FontSize = 20;
eval(sprintf('print -f%d -dpng MEG.png',f1))
eval(sprintf('print -f%d -dpng MEGcon.png',f2))

%% Analysis 2 (Comparison of frequency bands using variance/covariance)

clear
% Import data and define input cell array
load ('MEGPLANAR.mat'); % GRD
load ('y.mat')
Nband = numel(covariance)
V = cell(1,Nband+1);
V{Nband+1} = {[]};
for k=1:Nband
    V{k} = {covariance{k}};
    %V{k} = {variance{k}};
    if k==1
        V{Nband+1} = V{k};
    else
        V{Nband+1} = {V{Nband+1}{:},V{k}{1}};
    end
end
%V{end+1} = {V{1}{1},V{2}{1},V{3}{1},}; % All bands

% Classification step
rng('default') % For reproducibility
acc = mkl_class(V,y,'machine','easy',...
    'hyper',1,'CVratio',[0.8 0.2],...
    'Nrun',1000,'PCA_cut',0,'norm',1);

save FrqBnd_GrdCov acc
%save FrqBnd_GrdVar acc

% Plot resluts (Classification accuracy and Pos-hoc comparison)
titles = {'Delta','Theta','Alpha','Beta','lGamma','hGamma','All'};
pos_titles = {'Theta>Delta','Alpha>Theta','Beta>Alpha','lGamma>Beta','hGamma>lGamma','All>hGamma'};
%  define contrasts
c = [-1 1 0 0 0 0  0;
     0 -1 1 0 0 0  0;
     0 0 -1 1 0 0  0;
     0 0 0 -1 1 0  0;
     0 0 0 0 -1 1  0;
     0 0 0 0 0 -1  1
    ];
[f1,f2] = plot_results(titles,acc,pos_titles,c);
sgtitle('GRD COV')

eval(sprintf('print -f%d -dpng FBands.png',f1))
eval(sprintf('print -f%d -dpng FBandscon.png',f2))

%% Analysis 3 (MEG-MRI combination)
% Import data and define input cell array
clear

participants = spm_load(fullfile(wd,'participants-imputed.tsv'));
load ROIdata; load y; load MEGPLANAR; load cons;
cof = cons;
MRI = ROIdata;
MEG = covariance{6}; % Since lgamma does best numerically

V = {{cof},{MRI},{MEG},...
    {cof,MRI},{cof,MEG},{MRI,MEG},...
    {cof,MRI,MEG}};

% Classification step
rng('default') % For reproducibility
acc = mkl_class(V,y,'machine','easy',...
    'hyper',1,'CVratio',[0.8 0.2],...
    'Nrun',1000,'PCA_cut',0,'norm',1);

save CofsMRIMEG_GrdCov acc

% Plot resluts (Classification accuracy and Pos-hoc comparison)

% define titls in same order of the input cell V
titles = {'COF','MRI','MEG','MRI,COF','MEG,COF','MEG,MRI','MEG,MRI,COF'};
pos_titles = {'MRI>COF','MEG>COF','MEG>MRI',...
            'MRI,COF>COF','MEG,COF>COF','MEG,MRI>MRI',...
            'MEG,MRI,COF>MRI,COF'};
%  define contrasts
c = [-1 1 0  0 0 0  0; 
    -1 0 1   0 0 0  0; 
    0 -1 1   0 0 0  0;
    
    -1 0 0   1 0 0  0;
    -1 0 0   0 1 0  0; 
    0 -1 0   0 0 1  0; 
    
    0 0 0   -1 0 0  1];

[f1,f2] = plot_results(titles,acc,pos_titles,c)
sgtitle('GRD COV')

eval(sprintf('print -f%d -dpng CofMRIMEG.png',f1))
eval(sprintf('print -f%d -dpng CofMRIMEGMEGcon.png',f2))
