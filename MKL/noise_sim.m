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

load ('y.mat') % Classification labels, ie MCI vs CON
titles = {'Noise','Signal','Noise,Signal','Noise+Signal'};
pos_titles = {'Noise,Signal>Noise','Noise,Signal>Signal','Noise,Signal>Noise+Signal'};
%  define contrasts
c   = [-1 0 1 0;
       0 -1 1 0;
       0 0 1 -1];

Nk = 1000;
for n = [1 1000] % Number of noise element(s)
    
    rng('default') % for reproducibility

    accuracy = nan(Nk,size(c,2));
    for k = 1:Nk % Number of noise realisations
        
        noise = randn(numel(y), n);
        V = {{noise},{Signal},{Signal,noise},{[noise Signal]}};
        
        % Classification step
        rng('default') % For reproducibility
        acc = mkl_class(V,y,'machine','easy',...
            'hyper',1,'CVratio',[0.9 0.1],...
            'Nrun',1,'PCA_cut',0,'norm',1);
        
        accuracy(k,:) = mean(mean(acc,3),1);
    end
        
    % Plot resluts (Classification accuracy and Pos-hoc comparison)
    [f1] = plot_results(titles,accuracy,pos_titles,c)
    sgtitle(sprintf('n = %d', n))
    
end

%% Noise Simulation (new) - offered by Rik in the paper

clear titles y V x acc f
load ('y.mat') % Classification labels, ie MCI vs CON

% Define titles
titles = {'1 Noisy model (n=1)','1 Signal model (n=1)','1 Signal model (n=1) , 1000 Noisy models (n=1)','Noise (n=1000) + Signal (n=1) in single model'};
%pos_titles = {'Noise,Signal>Noise','Noise,Signal>Signal','Noise,Signal>Noise+Signal'};
%  define contrasts
%   c = [-1 0 1 0;
%       0 -1 1 0;
%       0 0 1 -1];

% Define input and output of classifier
V = cell(1,4);
rng(10) % For reproducibility
x=cell(1,1001);
x{1} = Signal;
for k = 2:numel(x)
    x{k} = randn(numel(y),1);
end
rng(10) % For reproducibility
V = {{randn(numel(y),1)},{Signal},x,{[randn(numel(y),1000) Signal]}};

% Classification step
rng('default') % For reproducibility
acc = mkl_class_ens(V,y,'machine','easy',...
    'hyper',1,'CVratio',[0.8 0.2],...
    'Nrun',1000,'PCA_cut',0,'norm',1);

save noise acc
% Plot resluts (Classification accuracy and Pos-hoc comparison)

f = plot_results(titles,acc);
eval(sprintf('print -f%d -dpng noise_sim_new.png',f))
