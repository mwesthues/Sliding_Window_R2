# Sliding_Window_R2

Initially, make this directory your home directory for running the analyses.
By default, this directory contains the subdirectory "/Data" and this 
README file. Please download the folder with R-scripts from GitHub and include 
them in the subdirectory "/Scripts" with the daughter directories "/Analysis"
and "/Functions". Further details are provided below.


# /Data
All raw data that are required for successfully running the scripts are stored
in "/Data".

# /Derived
All intermediate and final results will be stored in this folder.

# /Scripts
## Scripts/Analysis
Run the scripts in the following order:
(i) marker_preparation.R
(ii) r2v_sliding_windows.R
(iii) r2v_chr1.R
Please note that all scripts can be run in the terminal using "Rscript". 
Packages, which have not yet been installed by the user, will be automatically
installed and loaded. This requires the user to have specified a download 
mirror for CRAN packages.

## /Scripts/Functions
Scripts with functions that are called by scripts that are stored in 
"/Scripts/Analysis/".
(i) Marker_Functions.R
(ii) Misc_Functions.R
