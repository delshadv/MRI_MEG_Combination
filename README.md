# MRI_MEG_Combination
Combining MEG and MRI  to classify MCI

This repository contains all MATLAB codes including pre-processing, feature extraction and machine learning steps to reproduce results of paper ().
The dataset used for this project is "BioFIND"; the largest MEG data to study dementia which is freely available from DPUK (). For more info about data please see ()

It is necessary that you download and set up OSL toolbox from (https://ohba-analysis.github.io/osl-docs/) which consequently includes SPM toolbox (https://www.fil.ion.ucl.ac.uk/spm/). Depending on how you download OSL, you need to rename the OSL's root directory (after decompression) to "osl" and put it in your "current directory". This should be a folder called osl with at least the following contents:

osl/ . . .

├── MEG-ROI-nets/

├── osl-core/

├── parcellations/

├── spm12/ . . .

Then, clone this repository to your "current directory" with same level of osl (as below):

bwd/

├── osl/

├── MKL/ . . .

Note: You probably need MATLAB's "Parallel Computing" Toolbox to run the main codes within this repository faster.

"MKL" directory contains below scripts, functions, folder and files:

preproc_meg.m : MEG preprocessing steps  
preproc_mri.m : MRI preprocessing steps (DARTEL-VBM)  
participants-imputed.tsv : Participants info  
roi_extract.m : Extract ROI data from image volume(s)  
mkl_class.m : The function for Multi Kernel Learning for Binary Classification  
noise_sim.m : To reproduce supplementary results (Noise Simulation)  
main.m : Carries out permutation test for MCI vs Control classification problem to combine MEG - MRI  
HarvardOxford-combo-maxprob-thr25-2mm.nii : Harvard-Oxford atlas (116 ROIs)  
ctx_dist.csv :  Additional confound (mean distance from cortex to sensor)  
supplementary : This directory includes all complementary function and scripts needed to run above codes  

Work in progress

rik.henson@mrc-cbu.cam.ac.uk
