# MRI_MEG_Combination
Combining MEG and MRI  to classify MCI

This repository contains MATLAB codes including pre-processing, feature extraction and machine learning steps to reproduce results of paper (Late Combination shows that MEG adds to MRI in classifying MCI versus Controls published on Neuroimage : https://www.sciencedirect.com/science/article/pii/S1053811922001835?via%3Dihub).
The dataset used for this project is "BioFIND"; the largest MEG data to study dementia which is freely available from DPUK (https://portal.dementiasplatform.uk/Apply). For more info about data please see the dapaper on Neuroimage (https://www.sciencedirect.com/science/article/pii/S1053811922004633)

The codes are also available on DPUK once you applied for there. There is no need to download any toolbox through DPUK. However, if you want to analyze data out of DPUK portal please follow below:

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
mkl_ens_n.m : The function for late and intermediate combinations 
noise_sim.m : To reproduce supplementary results (Noise Simulation;first part (S1))   
supp2.m : To reproduce supplementary results (second part i.e. S2-S7)  
main.m : Carries out permutation test for MCI vs Control classification problem to combine MEG - MRI  
HarvardOxford-combo-maxprob-thr25-2mm.nii : Harvard-Oxford atlas (116 ROIs)  
ctx_dist.csv :  Additional confound (mean distance from cortex to sensor)  
supplementary : This directory includes all complementary functions and scripts, needed to run above codes  
derived : All csv files needed to reproduce results  

This code is for replicating experiment by paper titled "Late Combination shows that MEG adds to MRI in classifying MCI versus Control?"

Cite this paper if you are using part of code:

"Late Combination shows that MEG adds to MRI in classifying MCI versus Control" : https://www.sciencedirect.com/science/article/pii/S1053811922001835?via%3Dihub   

rik.henson@mrc-cbu.cam.ac.uk
