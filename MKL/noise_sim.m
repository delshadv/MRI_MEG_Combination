%% Supplementary results
% Combining MEG and MRI using MKL to classify MCI vs Control
% (BioFIND dataset)
%
% This script contains MEG preprocessing steps for reproducing results of the
% paper ""

% Henson R.N 2020, Vaghari D 2020

%% Define Paths ands variables

% Assumed you are currently in the directory including BioFIND data,
% OSL and MKL directories as described in readme.md
clear
%restoredefaultpath
cd MKL
addpath('supplementary'); 

participants = readtable('participants-imputed.tsv','FileType','text');
mri_num      = grp2idx(participants.sImaging);
%Signal = normalize(participants.Edu_years); % Education years
Signal = normalize(participants.MMSE); 

% Remove noisy-MRIs and non-MRI subjects
mri_num([23 197]) = 2;
Signal(mri_num==2,:) = [];

%% Noise simulation (part 1)
% Import data and define input cell array

load ('y.mat') % Classification labels, ie MCI vs CON
titles = {'Noise','Signal','Noise,Signal','Noise+Signal'};
%pos_titles = {'Noise,Signal>Noise','Noise,Signal>Signal','Noise,Signal>Noise+Signal'};
%  define contrasts
%   c = [-1 0 1 0;
%       0 -1 1 0;
%       0 0 1 -1];

for n = [1 10 100 1000 10000] %number of noise element(s)
    
    V = cell(1,4);
    rng('default') % For reproducibility
    noise = randn(numel(y), n);
    V = {{noise},{Signal},{Signal,noise},{[noise Signal]}};
    
    % Classification step
    rng('default') % For reproducibility
    acc = mkl_class(V,y,'machine','easy',...
        'hyper',1,'CVratio',[0.9 0.1],...
        'Nrun',1000,'PCA_cut',0,'norm',1);
    
    % Plot resluts (Classification accuracy and Pos-hoc comparison)

    [f1] = plot_results(titles,acc)
    sgtitle(sprintf('n = %d', n))
    
end

%% Noise simulation (part 2)
% Import data and define input cell array

n = 1; % number of noise
acc= []; accuracy = nan(1,1000);

for k=1:1000 % number of repetitions
   
    noise = randn(numel(y), n);
    V = {{noise}};
    
    % Classification step
    acc = mkl_class(V,y,'machine','easy',...
        'hyper',1,'CVratio',[0.9 0.1],...
        'Nrun',1,'PCA_cut',0,'norm',1);
    accuracy(k) = mean(acc,3);
end

figure('units','normalized','outerposition',[0 0 1 1]);
ax1 = subplot(3,1,1);histogram(accuracy); title('n=1 , r=1000');

n = 10; % number of noise
acc= []; accuracy = nan(1,100);

for k=1:100 % number of repetitions
   
    noise = randn(numel(y), n);
    V = {{noise}};
    
    % Classification step
    acc = mkl_class(V,y,'machine','easy',...
        'hyper',1,'CVratio',[0.9 0.1],...
        'Nrun',1,'PCA_cut',0,'norm',1);
    accuracy(2,k) = mean(acc,3);
end

ax2 = subplot(3,1,2);histogram(accuracy); title('n=10 , r=100');

n = 100; % number of noise
acc= []; accuracy = nan(1,10);

for k=1:10 % number of repetitions
   
    noise = randn(numel(y), n);
    V = {{noise}};
    
    % Classification step
    acc = mkl_class(V,y,'machine','easy',...
        'hyper',1,'CVratio',[0.9 0.1],...
        'Nrun',1,'PCA_cut',0,'norm',1);
    accuracy(k) = mean(acc,3);
end

ax3 = subplot(3,1,3);histogram(accuracy); title('n=100 , r=10');
linkaxes([ax1,ax2,ax3],'x')
 
