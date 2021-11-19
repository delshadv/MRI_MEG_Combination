function [bLcomI,bLcomL,acc1,acc2] = mkl_ens_nestedcv(V,y,varargin)
% [bLcomI,bLcomL,acc1,acc2] = mkl_ens(V,labels,machine,T,Nfold_i,Nfold_o,Nrun,PCA_cut,norm)
% Multi Kernel Learning for Binary Classification
%
% ens = flag for whether to do second stage of late combination of kernels
% (ensemble learning)
%
% acc1 : classification accuracy for "r" repetition with intermediate combination
% acc2 : classification accuracy for "r" repetition with late combination
% bLcomI : Hyperparameters for Intermediate combination (First level classifier)
% bLcomI : Hyperparameters for late combination (Second level classifier)
% y: targets (an array in which elements must be 1 for first class and 2
% for second class)
% V: inputs (a cell array, each element of cell can be features from a
% source of data e.g., MEG features, MRI features, etc and/or different
% combination of them)
% =================================================================================================================
% Optional args
% 'machine' : machine type, a string array ("pronto" for simpleMKL or
% "easy" for easyMKL ), default: "easy"

% 'T' : set threshold for feature selection based on MATLAB's rankfeature
%  function, default : None

% 'Nfold_o' : The Cross-Validation folds, default : 5 Folds CV

% 'Nfold_i' : The Cross-Validation folds, default : 10 Folds CV

% 'Nrun' : "r" repetition, default : 1000

% 'PCA_cut' : the percentage of the total variance explained by principal
% components, default : None (0: PCA off)

% 'norm' : 1 when normalization is needed, default : 1
% =================================================================================================================
% Example : mciConAccuracy = mkl_class(V,y,'machine','easy','PCA_run',95,'Nrun',10)
% =================================================================================================================
% classification functions were taken from the EasyMKL MATLAB implementation by Okba Bekhelifi.
% References:
% [1] Fabio Aiolli and Michele Donini
%      EasyMKL: a scalable multiple kernel learning algorithm
%      Paper @ http://www.math.unipd.it/~mdonini/publications.html
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
PCA_cut = 0;
feat_norm = 1;
%CVratio = [0.8 0.2];
Nfold_o = 5;
Nfold_i = 10;
Nrun = 10;
T = [];
ens = 1;

% Change default variables if they have been already defined by user
names = varargin(1:2:end);
values = varargin(2:2:end);

for k = 1:numel(names)
    %    names{k} = validatestring(names{k},{'machine','PCA_cut','norm','Nrun','Hyper','T','CVratio'}); % to ignore typos
    names{k} = validatestring(names{k},{'PCA_cut','Nfold_o','feat_norm','Nrun','T','Nfold_i','ens'}); % to ignore typos
    switch names{k}
        %         case 'machine'
        %             machine = values{k};
        
        case 'PCA_cut'
            PCA_cut = values{k};
            %        case 'CVratio'
            %            CVratio = values{k}; % It is better to import KFolds instead ; Done
        case 'Nfold_o'
            Nfold_o = values{k};
        case 'Nfold_i'
            Nfold_i = values{k};
        case 'feat_norm'
            feat_norm = values{k};
        case 'Nrun'
            Nrun = values{k};
        case 'T'
            T = values{k};
        case 'ens'
            ens = values{k};
    end
end

%% Balanced Cross-Validation

%Nfold = round(1/CVratio(2));
CVratio = [(Nfold_o-1)/Nfold_o 1-((Nfold_o-1)/Nfold_o)];
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

acc1 = cell(Nrun,1);
acc2 = cell(Nrun,1);
bLcomI = cell(Nrun,1);
bLcomL = cell(Nrun,1);

tracenorm = 1;

% test whether parpool  already open with numel(group etc...)
if ~isempty(gcp('nocreate'))
    try parpool(min(Nrun,12)); end
end

%tic
parfor rr = 1:Nrun
    
    rp = {};
    for c = 1:2
        rp{c} = randperm(N(c));
    end
    
    for aa=1:length(V)
        x = V{aa};
        
        for f = 1:Nfold_o
            
            % Spliting data from each dataset and create related kernel
            Ks_tr1 = []; Ks_ts1 = [];
            Ks_tr2 = []; Ks_ts2 = [];
            Ks_tr1C = []; Ks_ts1C = [];
            Ks_tr2C = []; Ks_ts2C = [];
            lambda1_f = []; %lambda2_f = [];
            Ks_tr1CC = {}; Ks_ts1CC = {};
            Ks_tr2CC = {}; Ks_ts2CC = {};
            ytrainCC = {} ;ytestCC = {} ;
            
            
            for qq = 1:length(x)
                
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
                
                % Best Hyper-parameter
                % Define and initialize variables
                gridC = [.01 .05 .1 .15 .2 .25 .3 .35 .4 .45 .5 .55 .6 .65 .7 .75 .8 .85 .9 .95 1];
                bestL1 = 0.01;
                LowestAcc1 = 0;
                
                for lambda1 = gridC
                    
                    CVratioC = [(Nfold_i-1)/Nfold_i 1-((Nfold_i-1)/Nfold_i)];
                    NC = [];
                    for c = 1:2
                        NC(c) = length(find(ytrain==c));
                        NtrainC(c)  = floor(CVratioC(1)*NC(c)); % was "ceil" - needs to be "floor"
                        NtestC(c) = NC(c) - NtrainC(c);
                    end
                    
                    % Ensure trained on equal number of each category
                    balNC = min(NtrainC); fldNC = floor(balNC*CVratioC(2));
                    for c = 1:2
                        if NtrainC(c)>balNC
                            NtestC(c) = NtestC(c)+(NtrainC(c)-balNC);
                            NtrainC(c) = balNC;
                        end
                    end
                    
                    
                    rpC = {};
                    for c = 1:2
                        rpC{c} = randperm(NC(c));
                    end
                    
                    L1 = ones(1, Nfold_i)*100;
                    %L2 = ones(1, Nfold_i)*100;
                    
                    for kC = 1:Nfold_i
                        
                        xtrainC = []; xtestC = []; ytrainC = []; ytestC = [];
                        
                        for c = 1:2
                            ii = setdiff([1:(balNC+fldNC)],[1:fldNC]+(kC-1)*fldNC);
                            ri = rpC{c}(ii);
                            rj = setdiff([1:NC(c)],ri);
                            
                            ii = find(ytrain==c);
                            %tr{f} = [tr{f}; ii(ri)]; te{f} = [te{f}; ii(rj)];
                            
                            xtrainC = [xtrainC; xtrain(ii(ri),:)];
                            xtestC  = [xtestC;  xtrain(ii(rj),:)];
                            ytrainC = [ytrainC; ytrain(ii(ri))];
                            ytestC  = [ytestC;  ytrain(ii(rj))];
                        end
                        
                        NtrainC = size(xtrainC,1); % = 2*balN
                        NtestC  = size(xtestC,1);
                        %Nfeats = size(xtrain,2);
                        
                        if feat_norm
                            xtrainC = normalize(xtrainC,'zscore');
                            xtestC  = normalize(xtestC,'zscore');
                        end
                        
                        % feature selection over training set
                        if ~isempty(T)
                            IDX1=[];
                            [IDX1, ~] = rankfeatures(xtrainC', ytrainC','CRITERION','roc');
                            xtrainC = xtrainC(:,IDX1(1:T));
                            xtestC  = xtestC(:,IDX1(1:T));
                        end
                        
                        if PCA_cut > 0
                            [coeff,scores,~,~,explained,mu] = pca(xtrainC);
                            sum_explained = 0;
                            idx = 0;
                            while sum_explained < PCA_cut
                                idx = idx + 1;
                                sum_explained = sum_explained + explained(idx);
                            end
                            xtrainC = scores(:,1:idx);
                            xtestC = (xtestC-mu)*coeff(:,1:idx);
                            
                        end
                        
                        ytestC(ytestC~=1)=-1;
                        ytrainC(ytrainC~=1)=-1;
                        
                        % L1-normalisation of kernels
                        tmp = xtrainC * xtrainC';
                        Ks_tr1C(:,:,qq) = reshape(normalize(tmp(:),'range'),NtrainC,NtrainC);
                        
                        tmp = xtestC * xtrainC';
                        Ks_ts1C(:,:,qq) = reshape(normalize(tmp(:),'range'),NtestC,NtrainC);
                        
                        if ens
                            easymkl_model = easymkl_train(Ks_tr1C(:,:,qq), ytrainC', lambda1, tracenorm);
                            [~, score1]   = easymkl_predict(easymkl_model, Ks_tr1C(:,:,qq));
                            [~, score2]   = easymkl_predict(easymkl_model, Ks_ts1C(:,:,qq));
                            
                            tmp = (score1 * score1');
                            Ks_tr2C(:,:,qq) = reshape(normalize(tmp(:),'range'),NtrainC,NtrainC);
                            
                            tmp = (score2 * score1');
                            Ks_ts2C(:,:,qq) = reshape(normalize(tmp(:),'range'),NtestC,NtrainC);
                        end
                        
                        if qq>=2
                            Ks_tr1CC{kC} = Ks_tr1C;
                            Ks_ts1CC{kC} = Ks_ts1C;
                            Ks_tr2CC{kC} = Ks_tr2C;
                            Ks_ts2CC{kC} = Ks_ts2C;
                            ytrainCC{kC} = ytrainC;
                            ytestCC{kC} = ytestC;
                        end

                        
                        easymkl_model  = easymkl_train(Ks_tr1C, ytrainC', lambda1, tracenorm);
                        [ts_pred, ~]   = easymkl_predict(easymkl_model, Ks_ts1C);
                        L1(kC) = ((sum(ts_pred==1 & ytestC==1)/sum(ytestC==1))+...
                            (sum(ts_pred==-1 & ytestC==-1)/sum(ytestC==-1)))/2*100;

                    end
                    
                    L1 = mean(L1);
                    if L1 > LowestAcc1
                        LowestAcc1 = L1 ;
                        bestL1 = lambda1 ;
                        %                         bestL2 = bestL2_tmp;
                    end
                    lambda1_f(qq,:) = bestL1;
                    %                    lambda2_f(qq,:) = bestL2;
                    
                end
                
                Ntrain = size(xtrain,1); % = 2*balN
                Ntest  = size(xtest,1);
                
                if feat_norm
                    xtrain = normalize(xtrain,'zscore');
                    xtest  = normalize(xtest,'zscore');
                end
                
                % feature selection over training set
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
                
                % L1-normalisation of kernels
                tmp = xtrain * xtrain';
                Ks_tr1(:,:,qq) = reshape(normalize(tmp(:),'range'),Ntrain,Ntrain);
                
                tmp = xtest * xtrain';
                Ks_ts1(:,:,qq) = reshape(normalize(tmp(:),'range'),Ntest,Ntrain);
                
                if ens
                    easymkl_model = easymkl_train(Ks_tr1(:,:,qq), ytrain', lambda1_f(qq,:), tracenorm);
                    [~, score1]   = easymkl_predict(easymkl_model, Ks_tr1(:,:,qq));
                    [~, score2]   = easymkl_predict(easymkl_model, Ks_ts1(:,:,qq));
                    
                    tmp = (score1 * score1');
                    Ks_tr2(:,:,qq) = reshape(normalize(tmp(:),'range'),Ntrain,Ntrain);
                    
                    tmp = (score2 * score1');
                    Ks_ts2(:,:,qq) = reshape(normalize(tmp(:),'range'),Ntest,Ntrain);
                end
                
            end
            
            if size(Ks_ts1,3)>1
                %% Tune lambda 2
                %gridC = [.01 .05 .1 .15 .2 .25 .3 .35 .4 .45 .5 .55 .6 .65 .7 .75 .8 .85 .9 .95 1];
                bestL2 = 0.01;
                LowestAcc2 = 0;
                LowestAcc3 = 0;
                
                for lambda2 = gridC
                    
                    L2 = ones(1, Nfold_i)*100;
                    L2I = ones(1, Nfold_i)*100;
                    
                    for kC = 1:Nfold_i
                        
                        ytrainC = ytrainCC{kC};
                        ytestC = ytestCC{kC};
                        Ks_tr1C = Ks_tr1CC{kC};
                        Ks_ts1C = Ks_ts1CC{kC};
                        
                        Best_easymkl_model  = easymkl_train(Ks_tr1C, ytrainC', lambda2, tracenorm);
                        [ts_pred, ~]   = easymkl_predict(Best_easymkl_model, Ks_ts1C);
                        L2I(kC) = ((sum(ts_pred==1 & ytestC==1)/sum(ytestC==1))+...
                            (sum(ts_pred==-1 & ytestC==-1)/sum(ytestC==-1)))/2*100;
                        
                        if ens
                            Ks_tr2C = Ks_tr2CC{kC};
                            Ks_ts2C = Ks_ts2CC{kC};
                            Best_easymkl_model  = easymkl_train(Ks_tr2C, ytrainC', lambda2, tracenorm);
                            [ts_pred, ~]   = easymkl_predict(Best_easymkl_model, Ks_ts2C);
                            L2(kC) =((sum(ts_pred==1 & ytestC==1)/sum(ytestC==1))+...
                                (sum(ts_pred==-1 & ytestC==-1)/sum(ytestC==-1)))/2*100;
                        end
                    end
                    L2I = mean(L2I);
                    if L2I > LowestAcc3
                        LowestAcc3 = L2 ;
                        bestL2I = lambda2 ;
                    end
                    
                    if ens
                        L2 = mean(L2);
                        if L2 > LowestAcc2
                            LowestAcc2 = L2 ;
                            bestL2 = lambda2 ;
                        end
                    end
                    
                end
            else
                bestL2I = bestL1;
                bestL2 = 1;
            end
            
            
            %%
            Best_easymkl_model  = easymkl_train(Ks_tr1, ytrain', bestL2I, tracenorm);
            [ts_pred, ~]   = easymkl_predict(Best_easymkl_model, Ks_ts1);
            acc1{rr}(aa,f) = ((sum(ts_pred==1 & ytest==1)/sum(ytest==1))+...
                (sum(ts_pred==-1 & ytest==-1)/sum(ytest==-1)))/2*100;
            bLcomI{rr}(aa,f) = bestL2I;
            
            if ens
                Best_easymkl_model  = easymkl_train(Ks_tr2, ytrain', bestL2, tracenorm);
                [ts_pred, ~]   = easymkl_predict(Best_easymkl_model, Ks_ts2);
                %acc2{rr}(aa,f) = (sum(ts_pred==ytest)/length(ytest))*100;
                acc2{rr}(aa,f) =((sum(ts_pred==1 & ytest==1)/sum(ytest==1))+...
                    (sum(ts_pred==-1 & ytest==-1)/sum(ytest==-1)))/2*100;
                bLcomL{rr}(aa,f) = bestL2;
            end
            
        end
    end
end

acc1 = permute(cell2mat(permute(acc1,[3,2,1])),[3,1,2]);
acc2 = permute(cell2mat(permute(acc2,[3,2,1])),[3,1,2]);
bLcomI = permute(cell2mat(permute(bLcomI,[3,2,1])),[3,1,2]);
bLcomL = permute(cell2mat(permute(bLcomL,[3,2,1])),[3,1,2]);

end

