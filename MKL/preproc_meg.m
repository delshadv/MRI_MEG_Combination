% Step 1
%% Pre-Processing of MEG data and defining confounds

% Combining MEG and MRI to classify MCI vs Control
% (BioFIND dataset)
%
% This script contains MEG preprocessing steps for reproducing results of the
% paper ""

% Henson R.N 2020, Vaghari D 2020

%% Define Paths ands variables

% Assumed you are currently in the directory including BioFIND data,
% OSL and MKL directories as described in readme.md

%restoredefaultpath
bwd = pwd;
wd  = fullfile(bwd,'MKL');
addpath(wd)

% Setup OSL
addpath(fullfile(bwd,'osl','osl-core'))
osl_startup
osl_check_installation

% BIDS and Processed directories
bidspth = fullfile(bwd,'BioFIND','MCIControls'); %BIDS Path
BIDS   = spm_BIDS(bidspth); % (If it does not work with OSL's SPM, so copy last version of spm_BIDS)
subs   = spm_BIDS(BIDS,'subjects', 'task', 'Rest');
nsub   = numel(subs);
subdir = cellfun(@(s) ['sub-' s], subs, 'UniformOutput',false);
procpth = fullfile(bidspth,'derivatives','meg_derivatives'); % If want maxfiltered files

% Define participant variables
participants = spm_load(fullfile(wd,'participants-imputed.tsv'));
group_num    = grp2idx(participants.group);
site_num     = grp2idx(participants.site);
sex_num      = grp2idx(participants.sex);
mri_num      = grp2idx(participants.sImaging);

% Remove noisy-MRIs and non-MRI subjects
mri_num([23 197]) = 2;
y = group_num(mri_num==1);  % Group labels
save(fullfile(bwd,'MKL','y'),'y')

%% Create Processed directory
% Please do all analysis in a separate directory from BIDS
% Here, we call it "Processed"

processed_pth = fullfile(bwd,'Processed');

if ~exist(processed_pth,'dir')
    
    mkdir('Processed');
    cd ('Processed')
    for s=1:nsub
        mkdir(sprintf('sub-Sub%04d',s))
    end
end

cd (processed_pth)

%% PreProcess- Part 1 (Convert, Downsample, Filter)

parfor sub = 1:nsub
    
    % Read event & json file to extract desirable length of MEG Recordings 
    tmp = spm_jsonread(fullfile(procpth,subdir{sub},'ses-meg1','meg',[subdir{sub} '_ses-meg1_task-Rest_proc-sss_meg.json']));
    event_file = spm_load(fullfile(bidspth,subdir{sub},'ses-meg1','meg',[subdir{sub} '_ses-meg1_task-Rest_events.tsv']));
    onset = (event_file.onset*tmp.SamplingFrequency)+1;
    
    % offset = onset + event_file.duration*tmp.SamplingFrequency;
    offset = (onset + 120 *tmp.SamplingFrequency)-1; % we put 120 seconds due to min length of raw data
    
    % Converting
    S = [];
    S.outfile = fullfile(processed_pth,subdir{sub},'spmeeg');
    S.dataset = fullfile(procpth,subdir{sub},'ses-meg1','meg',[subdir{sub} '_ses-meg1_task-Rest_proc-sss_meg.fif']);
    S.mode = 'epoched';
    S.channels = {'EOG', 'ECG', 'MEGMAG', 'MEGGRADPLANAR'}; % EEG was removed
    S.checkboundary = 0;
    S.trl = [onset offset 0];
    try
        S.conditionlabels = event_file.stim_type;
    catch
        S.conditionlabels = event_file.trial_type;
    end
    D = spm_eeg_convert(S);
    
    % Set channel types and bad channels
    S = [];
    S.D    = D;
    S.task = 'bidschantype';
    S.save = 1;
    S.filename = fullfile(bidspth,subdir{sub},'ses-meg1','meg',[subdir{sub} '_ses-meg1_task-Rest_channels.tsv']);
    D = spm_eeg_prep(S);
    D = chantype(D,indchantype(D,'MEGMAG'),'MEGMAG');
    D = chantype(D,indchantype(D,'MEGGRADPLANAR'),'MEGPLANAR');
    D.save
    
    % Downsampling the data
    S = [];
    S.D = D;
    S.method = 'resample';
    S.fsample_new = 500;
    D = spm_eeg_downsample(S);
    delete(S.D)
    
    % High-pass filter
    S = [];
    S.D = D;
    S.type = 'butterworth';
    S.band = 'high';
    S.freq = 0.5; % Cutoff frequency
    S.dir = 'twopass';
    S.order = 5;
    S.prefix = 'f';
    D = spm_eeg_filter(S);
    delete(S.D)
    
    
    % Low-pass filter
    S = [];
    S.D = D;
    S.type = 'butterworth';
    S.band = 'low';
    S.freq = 98; % Cutoff frequency
    S.dir = 'twopass';
    S.order = 5;
    S.prefix = 'f';
    D = spm_eeg_filter(S);
    delete(S.D)
    
end


%% PreProcess- Part 2 - Epoching, OSL Artifacts detection

parfor sub = 1:nsub
    
    infile = fullfile(processed_pth,subdir{sub},'effdspmeeg');
    D = spm_eeg_load(infile);
    
    EpochLength = 2 * D.fsample; % 2s in samples
    t = [1:EpochLength:(D.nsamples-EpochLength)]';
    nt = length(t);
    
    S = [];
    S.D = D;
    S.trl = [t t+EpochLength-1 zeros(nt,1)];
    S.conditionlabels = repmat({'EYES_CLOSED'},1,nt);
    S.prefix='e';
    S.bc = 0;
    D = spm_eeg_epochs(S);
    D.save;
    % OSL artifact detection
    D = osl_detect_artefacts(D,'modalities',unique(D.chantype(D.indchantype('MEGANY'))),'badchannels',false);
    D.save;
    
end

%% PreProcess- Part 3 - Despiking and calculate covariance matrix

freqbands = {[2 4],[4 8],[8 12],[12 30],[30 48],[52 86]};
modal = {'MEGPLANAR','MEGMAG'};

for k=1:numel(modal)
    
    covariance = cell(1,numel(freqbands));
    variance = cell(1,numel(freqbands));   
    
    for ii = 1:length(freqbands)
        
        Cov = cell(1,nsub); Var = Cov;
        
        parfor sub = 1:nsub
            
            infile = fullfile(processed_pth,subdir{sub},'effdspmeeg');
            D = spm_eeg_load(infile);
            
            % Remove bad badtrials
            chans = D.indchantype(modal{k},'GOOD'); % MEG :MAG or MEGPLANAR : GRD
            g = D(chans,:,:);
            g(:,:,D.badtrials)=[];
            
            % Generalise to Nchannels,Ntimes*Ntrials (good ones)
            y0 = reshape(g,size(g,1),D.nsamples*size(g,3));
            
            % Filter to desired freq band
            y1 = ft_preproc_bandpassfilter(y0, D.fsample, freqbands{ii}, 4, 'but');
            
            % Despiking
            y = filloutliers(y1','clip','median','ThresholdFactor',3);
            
            % Calculate Covariance Matrix
            cm = cov(y);
            
            Cov{sub} = cm(find(triu(cm,1)))';
            Var{sub} = diag(cm)';
            
        end
        
        Cov = cat(1,Cov{:});
        Var = cat(1,Var{:});
        covariance{ii} = Cov(mri_num==1,:);
        variance{ii} = Var(mri_num==1,:);
    end
    
    save(sprintf('%s',modal{k}), 'covariance','variance');
    
end

%% Define and save Confounds

ctx_dist = csvread(fullfile(bwd,'MKL','ctx_dist.csv'));
cofs = [normalize(participants.age) normalize(participants.Move1) ...
    normalize(participants.Move2) normalize(participants.Recording_time) ...
    normalize(participants.sex) normalize(participants.site) ...
    normalize(ctx_dist(:,2).^2) normalize(participants.Edu_years)];
cons(mri_num==2,:) = [];
save(fullfile(bwd,'MKL','cofs'),'cofs')



