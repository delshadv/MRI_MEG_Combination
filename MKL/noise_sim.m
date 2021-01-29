%% Supplementary results
% Combining MEG and MRI using MKL to classify MCI vs Control
% (BioFIND dataset)
%
% This script contains MEG preprocessing steps for reproducing results of the
% paper ""

% Henson R.N 2020, Vaghari D 2020

%% Define Paths ands variables

% Assumed you are currently in the directory including BioFIND data,
cd MKL
addpath('supplementary'); 
 
clear

participants = readtable('participants-imputed.tsv','FileType','text');
% Remove noisy-MRIs and non-MRI subjects
mri_num      = grp2idx(participants.sImaging);
mri_num([23 197]) = 2;
participants(mri_num==2,:) = [];

X = [participants.Edu_years participants.MMSE participants.age participants.sex]; 

y = csvread('derived/labels.csv'); % (better name for file?) What about labels?

No  = size(y,1) % Number of observations
Nsf = size(X,2) % Number of signal features
Nnf = 100       % Number of noise features

Np = 100;      % Number of  noise randomisations (not permutations)

% test whether parpool  already open with numel(group etc...)
if isempty(gcp('nocreate')) 
    try parpool(min(Np,12)); end
end

rng('default') % for reproducibility

Signal = X; 

accuracy1 = cell(1,Np); % For intermediate combination
accuracy2 = cell(1,Np); % For late combination
    
parfor p = 1:Np % Number of noise realisations
    
    Noise = randn(No, Nnf);
    
    V = {};
    V{1} = {Noise};
    V{2} = {Signal(:,1)};
    V{3} = {[Signal(:,1) Noise]}; % early combination (concatenation)
    % V{4} = {Signal(:,1), Noise};  % two separate kernels
    V{4} = mat2cell([Signal(:,1) Noise],No,ones(1,1+Nnf));  % all separate kernels
    V{5} = {Signal(:,1), Signal(:,2:4)};

   % Intermediate better when C1=1 than C1=0.1; Late better when C1=0.1 
   [acc1,~] = mkl_ens(V,y,'Hyper1',1,'CVratio',[0.8 0.2],'Nrun',1,'PCA_cut',0,'feat_norm',1,'ens',0);
   [~,acc2] = mkl_ens(V,y,'Hyper1',0.1,'Hyper2',1,'CVratio',[0.8 0.2],'Nrun',1,'PCA_cut',0,'feat_norm',1,'ens',1);
%    [acc1,acc2] = mkl_ens(V,y,'Hyper1',0.1,'Hyper2',1,'CVratio',[0.8 0.2],'Nrun',1,'PCA_cut',0,'feat_norm',1,'ens',1);
    
    accuracy1{p} = mean(mean(acc1,3),1);
    accuracy2{p} = mean(mean(acc2,3),1);
    fprintf('p=%d, acc1=%s, acc2=%s\n',p,mat2str(round(accuracy1{p})),mat2str(round(accuracy2{p})))
end

accuracy1 = cat(1,accuracy1{:});
mean(accuracy1)
accuracy2 = cat(1,accuracy2{:});
mean(accuracy2)

% Plot resluts (Classification accuracy and Pos-hoc comparison)

titles = {'Noise','Signal1','Signal1+Noise','Signal1,Noise','Signal1,Signal2-4'};
%  define contrasts
c   = [-1 1 0 0 0;
       0 1 0 -1 0;
       0 0 -1 1 0];
con_titles = {'Signal1>Noise','Signal1>Signal1,Noise','Signal1,Noise>Signal1+Noise'};

% [f1] = plot_results(titles,accuracy1,con_titles,c)
% sgtitle(sprintf('Intermediate, Nf = %d', Nnf))
% 
% [f1] = plot_results(titles,accuracy2,con_titles,c)
% sgtitle(sprintf('Late, Nf = %d', Nnf))

[f1] = plot_results(cat(2,titles,titles),[accuracy1 accuracy2],cat(2,titles,titles),[-eye(5) eye(5)])
sgtitle(sprintf('Intermediate+Late, Nf = %d', Nnf))


return

% Or show same effect on independent dataset...
% % Create 2 predictor example from Matlab's example data
% fish = load('/hpc-software/matlab/r2019a/toolbox/stats/statsdemos/fisheriris.mat')
% inds = ~strcmp(fish.species,'setosa');
% X = fish.meas(inds,:); 
% y = fish.species(inds);
% y = grp2idx(y);
% save('supplementary/fish_data','X','y')

%load('supplementary/fish_data')

