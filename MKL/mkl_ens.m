function [acc1,acc2] = mkl_ens(V,y,varargin)
% [acc1,acc2] = mkl_ens(V,labels,machine,hyper1,hyper2,T,CVratio,Nrun,PCA_cut,norm)
% Multi Kernel Learning for Binary Classification 
%
% ens = flag for whether to do second stage of late combination of kernels
% (ensemble learning)
%
% acc1 : classification accuracy for "r" repetition with intermediate combination
% acc2 : classification accuracy for "r" repetition with late combination
%
% y: targets (an array in which elements must be 1 for first class and 2
% for second class)
% V: inputs (a cell array, each element of cell can be features from a 
% source of data e.g., MEG features, MRI features, etc and/or different 
% combination of them)
% =================================================================================================================
% Optional args
% 'machine' : machine type, a string array ("pronto" for simpleMKL or
% "easy" for easyMKL ), default: "easy"

% 'Hyper1' : MKL lambda hyperparameter - penalty or regularization term, default: 0.1 for feature combination

% 'Hyper2' : MKL lambda hyperparameter - penalty or regularization term, default: 1 for decision combination

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
lambda1 = 0.1;
lambda2 = 1;
PCA_cut = 0;
feat_norm = 1;
CVratio = [0.8 0.2];
Nrun = 1000;
T = [];

% Change default variables if they have been already defined by user
names = varargin(1:2:end);
values = varargin(2:2:end);

for k = 1:numel(names)
%    names{k} = validatestring(names{k},{'machine','PCA_cut','norm','Nrun','Hyper','T','CVratio'}); % to ignore typos
    names{k} = validatestring(names{k},{'PCA_cut','feat_norm','Nrun','Hyper1','Hyper2','T','CVratio','ens'}); % to ignore typos
    switch names{k}
%         case 'machine'
%             machine = values{k};
        case 'Hyper1'
            lambda1 = values{k};
        case 'Hyper2'
            lambda2 = values{k};
        case 'PCA_cut'
            PCA_cut = values{k};
        case 'CVratio'
            CVratio = values{k}; % It is better to import KFolds instead
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

acc1 = cell(Nrun,1);
acc2 = cell(Nrun,1);

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
        
        for f = 1:Nfold
            
            % Spliting data from each dataset and create relted kernel
            Ks_tr1 = []; Ks_ts1 = [];
            Ks_tr2 = []; Ks_ts2 = [];
            
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
                 
                Ntrain = size(xtrain,1); % = 2*balN
                Ntest  = size(xtest,1);
                %Nfeats = size(xtrain,2);
                
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
                    easymkl_model = easymkl_train(Ks_tr1(:,:,qq), ytrain', lambda1, tracenorm);  
                    [~, score1]   = easymkl_predict(easymkl_model, Ks_tr1(:,:,qq));
                    [~, score2]   = easymkl_predict(easymkl_model, Ks_ts1(:,:,qq));
                    
                    tmp = (score1 * score1');
                    Ks_tr2(:,:,qq) = reshape(normalize(tmp(:),'range'),Ntrain,Ntrain); 
                    
                    tmp = (score2 * score1');
                    Ks_ts2(:,:,qq) = reshape(normalize(tmp(:),'range'),Ntest,Ntrain);    
                end
            end
             
            easymkl_model  = easymkl_train(Ks_tr1, ytrain', lambda1, tracenorm);
            [ts_pred, ~]   = easymkl_predict(easymkl_model, Ks_ts1);
            %acc1{rr}(aa,f) = (sum(ts_pred==ytest)/length(ytest))*100;
            acc1{rr}(aa,f) = ((sum(ts_pred==1 & ytest==1)/sum(ytest==1))+...
                (sum(ts_pred==-1 & ytest==-1)/sum(ytest==-1)))/2*100;
            if ens
                easymkl_model  = easymkl_train(Ks_tr2, ytrain', lambda2, tracenorm);
                [ts_pred, ~]   = easymkl_predict(easymkl_model, Ks_ts2);
                %acc2{rr}(aa,f) = (sum(ts_pred==ytest)/length(ytest))*100;
                acc2{rr}(aa,f) =((sum(ts_pred==1 & ytest==1)/sum(ytest==1))+...
                    (sum(ts_pred==-1 & ytest==-1)/sum(ytest==-1)))/2*100;
            end     
        end
    end
end
%toc;

% Better way of doing below!? Yes, using permute
%tmp1 = []; tmp2 = [];
%for rr = 1:Nrun
%    tmp1(rr,:,:) = acc1{rr};
%    tmp2(rr,:,:) = acc2{rr};
%end
%acc1 = tmp1; acc2 = tmp2;
acc1 = permute(cell2mat(permute(acc1,[3,2,1])),[3,1,2]);
acc2 = permute(cell2mat(permute(acc2,[3,2,1])),[3,1,2]);

end

