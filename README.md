# Github-MKL
Combining MEG and MRI using Multiple Kernel Learning to classify MCI

This repository contains all MATLAB scripts (and functions) including pre-processing, feature extraction and machine learning steps to reproduce results of paper ().
Dataset (BioFIND) used for this projects are freely available from DPUK in (). For more info about data please see ()

Before downloading the repository: It is necessary that you first download and set up OSL toolbox from (https://ohba-analysis.github.io/osl-docs/) which consequently includes SPM toolbox (https://www.fil.ion.ucl.ac.uk/spm/). Depending on how you download OSL, you need to rename the OSL's root directory (after decompression) to "osl" and put it in your current directory. This should be a folder called osl with at least the following contents:

osl/ . . .

├── MEG-ROI-nets/

├── osl-core/

├── parcellations/

├── spm12/ . . .

Then, clone this repository to the root directory (same level with osl) and download all of this repository within code folder. So we will have a directory structure in your current directory like this:

bwd/

├── osl/

├── MKL/ . . .

Note: It is necessary to have MATLAB's "Digital Signal Processing" and "Statistics and Machine Learning" and "Parallel Computing " Toolboxes installed within the MATLAB version you use to run this repository successfully.

MKL repository content:

   

Work in progress

rik.henson@mrc-cbu.cam.ac.uk
