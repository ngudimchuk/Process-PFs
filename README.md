# Process-PFs
Matlab 2017b scripts for analysis of curled protofilaments (PFs), traced in IMOD software [refs: McIntosh et al., JCB 2018; Gudimchuk et al., Nature Commun., 2020]

Open Application_PF_processing.mlapp
- specify the name of your dataset
- specify a folder where your files with PF tracings are located (.dat in IMOD-generated format)
- specify a folder to store the output files with analysis results (by default a new folder "Results" will be created in the folder with the raw data points)
- check/uncheck box to proccess 3d protofilaments (yes) or 2d protofilaments (no)
- specify the pixel size
- check/uncheck boxes to view individual PFs before/after processing (press any button to continue or "close all figures" to stop)
- press "Start Processing" button
- (optionally): to exclude outlier protofilament tracings from analysis, create 'excludePFs.txt' file in the folder containing your .dat files. The txt file should contain a vector of numbers for the protofilaments to be excluded, e.g. 1 20 21 30

Output txt files:

- ang_vs_tip - dependence of average PF curvature on the distance from PF tip {distance(nm); curvature (deg/dimer); SD of curvature (deg/dimer)} 
- ang_vs_tip_fit - weighted linear fit of dependence of average PF curvature on the distance from PF tip {distance(nm)}
- all_angles - individual values of all angles between consecutive pairs of PF segments (deg/dimer) 
- Ang_ALL_mean_med_SD_N - average PF curvature (deg/dimer); median PF curvature (deg/dimer); SD of PF curvature (deg/dimer); N of angles 
- hist_ang_all - histogram of PF curvatures
- all_PF_lengths - individual values of PF lengths (nm)
- hist_PF_length - histogram of PF lengths
- MeanL_StdL_numL - mean PF length, SD pf PF length, N of PFs
- PL_SDPL - dynamic persistence length, SD of dynamic the persistence length
- dist_logcos_errlogcos - distance along the PF (nm) vs. the logarithm of the averaged cosine of the deflection of the protofilaments from their average shape vs. error of the measurement of the log(cos<teta>) vs. the ordinate of the weighted linear fit



Output jpeg/png files:
- angles - histogram of PF curvatures;
- histlegth - histogram of PF lengths;
- plotPFs - all PF traces, smoothed
- PL - logarithm of the averaged cosine of the deflection of the protofilaments from their average shape (blue) and a weighted linear fit (red)
