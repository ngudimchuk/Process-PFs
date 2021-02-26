Matlab 2017b scripts for analysis of curled protofilaments, traced in IMOD software [refs: McIntosh et al., JCB 2018; Gudimchuk et al., Nature Commun., 2020; XXX]

Open Application_PF_processing.mlapp
%specify the name of your dataset
%specify a folder where your files with PF tracings are located (.dat in iMod-generated format)
%specify a folder to store the output files with analysis results (by default a new folder "Results" is made in the folder with raw data points)
%specify the pixel size
%check/uncheck boxes to view individual PFs before/after processing (press any button to continue or "close all figures" to stop)
%press "Start Processing" button

Output txt files:

ang_vs_tip - dependence of average PF curvature on the distance from PF tip {distance(nm); curvature (deg/dimer); SD of curvature (deg/dimer)} 
ang_vs_tip_fit - weighted linear fit of dependence of average PF curvature on the distance from PF tip {distance(nm)} 
Ang_ALL_mean_med_SD_N - average PF curvature (deg/dimer); median PF curvature (deg/dimer); SD of PF curvature (deg/dimer); N of angles 
hist_ang_all - histogram of PF curvatures
hist_PF_length - histogram of PF lengths
MeanL_StdL_numL - mean PF length, SD pf PF length, N of PFs


Output jpeg/png files:
angles - histogram of PF curvatures
histlegth - histogram of PF lengths
plotPFs - all PF traces, smoothed


