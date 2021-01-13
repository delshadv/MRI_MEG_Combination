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

%%Signal = normalize(participants.MMSE); 
%Signal = normalize([participants.MMSE participants.Edu_years]); 
%Nsf = size(Signal,2)

Nsf = 100;
Signal = repmat(participants.MMSE,1,Nsf)+randn(size(participants.MMSE,1),Nsf);

% Remove noisy-MRIs and non-MRI subjects
mri_num([23 197]) = 2;
Signal(mri_num==2,:) = [];


load ('y.mat') % Classification labels, ie MCI vs CON
titles = {'Noise','Signal1','Signal1,Signal2','Signal1,Noise','Signal1+Noise'};
%  define contrasts
c   = [0 -1 1 0 0;
       0 1 0 -1 0;
       0 0 0 1 -1];
con_titles = {'Signal1,Signal2>Signal1','Signal1>Signal1,Noise','Signal1,Noise>Signal1+Noise'};

Np = 100;
parpool(min(Np,12));

for Nnf = 100 %[1 100] % Number of noise feature(s)
    
    rng('default') % for reproducibility

    accuracy1 = cell(1,Np); accuracy2 = accuracy1;
 
    parfor p = 1:Np % Number of noise realisations
        
        Noise = randn(numel(y), Nnf);
        sepS  = mat2cell(Signal,numel(y),ones(1,Nsf));
        sepSN = mat2cell([Signal(:,1) Noise],numel(y),ones(1,1+Nnf));
        V = {{Noise},{Signal(:,1)},sepS,sepSN,{[Signal(:,1) Noise]}};
        
        [acc1,acc2] = mkl_ens(V,y,'Hyper1',0.1,'Hyper2',0.1,'CVratio',[0.8 0.2],'Nrun',1,'PCA_cut',0,'feat_norm',1,'ens',1);

        %% faster, if late = intermediate when only 1 kernel...?
        %         V = {{noise},{Signal},sepks,{[noise Signal]}};
        %         [acc1,~] = mkl_ens(V,y,'hyper',1,'CVratio',[0.8 0.2],'Nrun',1,'PCA_cut',0,'feat_norm',1,'ens',0);
%         V = {sepks};
%         [~,acc2] = mkl_ens(V,y,'hyper',1,'CVratio',[0.8 0.2],'Nrun',1,'PCA_cut',0,'feat_norm',1,'ens',1);
        
        accuracy1{p} = mean(mean(acc1,3),1);
        accuracy2{p} = mean(mean(acc2,3),1);
        fprintf('p=%d, acc1=%s, acc2=%s\n',p,mat2str(round(accuracy1{p})),mat2str(round(accuracy2{p})))
    end
        
    % Plot resluts (Classification accuracy and Pos-hoc comparison)
     accuracy1 = cat(1,accuracy1{:});
%     [f1] = plot_results(titles,accuracy1,con_titles,c)
%     sgtitle(sprintf('Intermediate, Nf = %d', Nnf))
%   
     accuracy2 = cat(1,accuracy2{:});
%     [f1] = plot_results(titles,accuracy2,con_titles,c)
%     sgtitle(sprintf('Late, Nf = %d', Nnf))

    [f1] = plot_results(cat(2,titles,titles),[accuracy1 accuracy2],cat(2,titles,titles),[-eye(5) eye(5)])
    sgtitle(sprintf('Intermediate+Late, Nf = %d', Nnf))

end


