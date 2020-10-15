
% Step 2
%% Pre-Processing of MRI data

% Combining MEG and MRI using MKL to classify MCI vs Control
% (BioFIND dataset)
%
% This script contains MRI preprocessing steps for reproducing results of the
% paper ""

% Henson R.N 2020, Vaghari D 2020

%% Define Paths ands variables

% Assumed you are currently in the directory including "BioFIND" data,
% "OSL" and "MKL directories as described in readme.md (current folder)

%restoredefaultpath
bwd = pwd;
wd  = fullfile(bwd,'MKL');
addpath(fullfile(bwd,'osl','osl-core'))
addpath(fullfile(bwd,'MKL'))
addpath(fullfile(bwd,'MKL','supplementary'))
osl_startup
osl_check_installation

bidspth = fullfile(bwd,'BioFIND','MCIControls'); %BIDS Path
BIDS   = spm_BIDS(bidspth); % (If it does not work with OSL's SPM, so copy last version of spm_BIDS)
subs   = spm_BIDS(BIDS,'subjects', 'task', 'Rest');
nsub   = numel(subs);
subdir = cellfun(@(s) ['sub-' s], subs, 'UniformOutput',false);

spm('defaults', 'FMRI');
spm_jobman('initcfg')

%% Create Processed_mri directory
% Please do all analysis in a separate directory from BIDS
% Here, we call it "Processed_mri"
processed_pth= fullfile(bwd,'Processed_mri'); %for testing: /imaging/dv01/defacing_project/MRI_defaced
if ~exist(processed_pth,'dir')
    mkdir(processed_th);  
end
cd (processed_pth)

for sub = 1:nsub
    if sub == 23 || sub == 127
        'Noisy MRI'
    else
        try
            mrfile = fullfile(bwd,'BioFIND','MCIControls',subdir{sub},'ses-meg1','anat',[subdir{sub} '_ses-meg1_T1w.nii.gz']);
            newfile = fullfile(processed_pth,[subdir{sub} '_ses-meg1_T1w.nii.gz']);
            copyfile(mrfile,newfile);
            eval(sprintf('gunzip %s',newfile));
            delete(newfile)           
        catch
            warning('No MRI')
        end
    end
end

P = spm_select('FPList',processed_pth,'^sub.*nii');

%% First job to segment each T1
jobfile = {fullfile(bwd,'MKL','supplementary','batch_seg_job.m')};

parfor s=1:size(P,1)
    inputs = {};
    inputs{1} = cellstr(P(s,:));
    spm_jobman('run', jobfile, inputs{:});
end

%% Second job combines images for DARTEL 

jobfile = {fullfile(bwd,'MKL','supplementary','batch_dartel_job.m')};

inputs = {};

P1 = spm_select('FPList',processed_pth,'^rc1.*nii');
inputs{1} = cellstr(P1);

P2 = spm_select('FPList',processed_pth,'^rc2.*nii');
inputs{2} = cellstr(P2); % assumes same order as above!

spm_jobman('run', jobfile, inputs{:});

%% Third job writes images to MNI space (can be combined with above if don't want parallel)

jobfile = {fullfile(bwd,'MKL','supplementary','batch_mni_job.m')};

inputs = {};

P1 = spm_select('FPList',processed_pth,'^Template.*nii');

P2 = spm_select('FPList',processed_pth,'^u_.*nii');

P3 = spm_select('FPList',processed_pth,'^c1s.*nii'); 
  
parfor s=1:size(P3,1)
    inputs = {};
    inputs{1} = cellstr(P1(end,:));
    inputs{2} = cellstr(P2(s,:));
    inputs{3} = cellstr(P3(s,:));
    spm_jobman('run', jobfile, inputs{:});    
end

%% Calculate ROI data and save it
S = [];
S.ROIfiles{1} = fullfile(bwd,'MKL','HarvardOxford-combo-maxprob-thr25-2mm.nii'); % 116 ROIs
P = spm_select('FPList',processed_pth,'^swc1sub.*nii');

for s=1:size(P,1)
    S.Datafiles{1}{s} = P(s,:);
end

ROI = roi_extract(S);

ROIdata = [];
for r = 1:length(ROI)
    ROIdata(:,r) = ROI(r).mean;
    %num_vox(:,r) = ROI(r).numvox;
end

ROIdata(:,111:end) = [];

save(fullfile(bwd,'MKL','ROIdata'),'ROIdata')


