
function [acc] = mkl_class_ens(V,y,varargin)
% [acc] = MKL_class(V,y,machine,C,T,CVratio,Nrun,PCA_cut,norm)
% Multi Kernel Learning for Binary Classification (it can be used for 
%regression models after a few modifications)
%
% acc : classification accuracy for "r" repetition
% y: targets (an array in which elements must be 1 for first class and 2
% for second class)
% V: inputs (a cell array, each element of cell can be features from a 
% source of data e.g., MEG features, MRI features, etc and/or different 
% combination of them)
% =================================================================================================================
% Optional args
% 'machine' : machine type, a string array ("pronto" for simpleMKL or
% "easy" for easyMKL ), default: "easy"

% 'Hyper' : MKL C hyperparameter - penalty or regularization term, default : 1

% 'T' : set threshold for feature selection based on MATLAB's rankfeature 
%  function, default : None

% 'CVratio' : The Cross-Validation ratio, default : [0.8 0.2] ==> 5 Folds CV

% 'Nrun' : "r" repetition, default : 1000

% 'PCA_cut' : the percentage of the total variance explained by principal
% components, default : None (0: PCA off)

% 'norm' : 1 when normalization is needed, default : 1
% =================================================================================================================
% Example : mciConAccuracy = mkl_class(V,y,'machine','easy','PCA_run',95,'Nrun',10)
% =================================================================================================================
% classification functions were taken from the PRONTO toolbox and EasyMKL MATLAB implementation by Okba Bekhelifi.
% by Rik Henson, Delshad Vaghari Sep,2020
% =================================================================================================================
%% Interface
if nargin < 2
    error('Not enough input arguments')
end

if mod(nargin,2)==1
    error('There must be an even number (pairwise) of inputs.')
end

% Define default variables
C = 1;
machine = "easy";
PCA_cut = 0;
norm = 1;
CVratio = [0.8 0.2];
Nrun = 1000;
T = [];

% Change default variables if they have been already defined by user
names = varargin(1:2:end);
values = varargin(2:2:end);

for k = 1:numel(names)
    names{k} = validatestring(names{k},{'machine','PCA_cut','norm','Nrun','Hyper','T','CVratio'}); % to ignore typos
    switch names{k}
        case 'machine'
            machine = values{k};
        case 'Hyper'
            C = values{k};
        case 'PCA_cut'
            PCA_cut = values{k};
        case 'CVratio'
            CVratio = values{k}; % It is better to import KFolds instead
        case 'norm'
            norm = values{k};
        case 'Nrun'
            Nrun = values{k};
        case 'T'
            T = values{k};
    end
end

%% Balanced Cross-Validation

tic;
Nfold = round(1/CVratio(2));
%rbf = @(X,Y) exp(-0.1 .* pdist2(X,Y, 'euclidean').^2); % RBF kernel

for c = 1:2
    N(c) = length(find(y==c));
    Ntrain(c)  = floor(CVratio(1)*N(c)); % was "ceil" - needs to be "floor"
    Ntest(c) = N(c) - Ntrain(c);
end

% Ensure trained on equal number of each category
balN = min(Ntrain); fldN = floor(balN*CVratio(2));
for c = 1:2
    if Ntrain(c)>balN
        Ntest(c) = Ntest(c)+(Ntrain(c)-balN);
        Ntrain(c) = balN;
    end
end

acc = cell(Nrun,1);

parfor rr = 1:Nrun
    
    rp = {};
    for c = 1:2
        rp{c} = randperm(N(c));
    end
    
    for aa=1:length(V)
        x = V{aa};
        
        for f = 1:Nfold
            
            % Spliting data from each dataset and create relted kernel
            d = struct();Ks_tr=[];Ks_ts=[];
            
            for qq=1:length(x)
                
                X = x{qq};
                xtrain = []; xtest = []; ytrain = []; ytest = [];
                for c = 1:2
                    ii = setdiff([1:(balN+fldN)],[1:fldN]+(f-1)*fldN);
                    ri = rp{c}(ii);
                    rj = setdiff([1:N(c)],ri);
                    
                    ii = find(y==c);
                    %tr{f} = [tr{f}; ii(ri)]; te{f} = [te{f}; ii(rj)];
                    
                    xtrain = [xtrain; X(ii(ri),:)];
                    xtest  = [xtest;  X(ii(rj),:)];
                    ytrain = [ytrain; y(ii(ri))];
                    ytest  = [ytest;  y(ii(rj))];
                end
                
                d.tr_targets = ytrain;
                d.te_targets = ytest;
                
                if (norm)
                    [xtrain,xtest]=normalizemeanstd(xtrain,xtest);
                end
                %feature selection over each individual dataset
                
                if ~isempty(T)
                    IDX1=[];
                    [IDX1, ~] = rankfeatures(xtrain', ytrain','CRITERION','roc');
                    xtrain = xtrain(:,IDX1(1:T));
                    xtest  = xtest(:,IDX1(1:T));
                end
                if PCA_cut > 0
                    [coeff,scores,~,~,explained,mu] = pca(xtrain);
                    sum_explained = 0;
                    idx = 0;
                    while sum_explained < PCA_cut
                        idx = idx + 1;
                        sum_explained = sum_explained + explained(idx);
                    end
                    xtrain = scores(:,1:idx);
                    xtest = (xtest-mu)*coeff(:,1:idx);
                    
                end
                ytest(ytest~=1)=-1;
                ytrain(ytrain~=1)=-1;
                
                tmp = (xtrain * xtrain');
                normTmp = tmp- min(tmp(:));
                d.train1 = normTmp ./ max(normTmp(:)); % *
                
                tmp = (xtest * xtrain');
                normTmp = tmp- min(tmp(:));
                d.test1 = normTmp ./ max(normTmp(:)); % *
                
                Ks_tr1 = d.train1;
                Ks_ts1 = d.test1
                tracenorm = 1;
                easymkl_model = easymkl_train(Ks_tr1, ytrain', 0.1, tracenorm);
                [~, score1] = easymkl_predict(easymkl_model, Ks_tr1);
                [~, score2] = easymkl_predict(easymkl_model, Ks_ts1);


                tmp = (score1 * score1');
                normTmp = tmp- min(tmp(:));
                d.train{qq}  = normTmp ./ max(normTmp(:)); % *
                tmp = (score2 * score1');
                normTmp = tmp- min(tmp(:));
                d.test{qq}  = normTmp ./ max(normTmp(:)); % *

                Ks_tr(:,:,qq) = d.train{qq};
                Ks_ts(:,:,qq) = d.test{qq};
            end
            
            
%% Classification Model
            
            switch(machine)
                
                case 'easy'
                    ytest(ytest~=1)=-1;
                    ytrain(ytrain~=1)=-1;
                    
                    tracenorm = 1;
                    easymkl_model = easymkl_train(Ks_tr, ytrain', C, tracenorm);
                    [ts_pred, ~] = easymkl_predict(easymkl_model, Ks_ts);
                    acc{rr}(aa,f) = (sum(ts_pred==ytest)/length(ytest))*100;
                    
                    
                case 'pronto'
                    opti = prt_machine_sMKL_cla(d,C);
                    acc{rr}(aa,f) = sum((opti.predictions==d.te_targets)/length(ytest))*100;
            end
        end
        
    end
    
end

toc;

tmp = [];
for rr = 1:Nrun
    tmp(rr,:,:) = acc{rr};
end
acc = tmp;
clear tmp
end

