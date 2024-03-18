# Process-PFs

MATLAB 2017b scripts for the analysis of curled protofilaments (PFs), traced in IMOD software [refs: McIntosh et al., JCB 2018; Gudimchuk et al., Nature Commun., 2020].

Open Application_PF_processing.mlapp:

- Specify the name of your dataset.
- Specify the folder where your files with PF tracings are located (.dat in IMOD-generated format).
- Specify a folder to store the output files with analysis results (if not specified, a new folder named "Results" will be created in the folder with the raw data points).
- Check/uncheck the box to process 3D protofilaments (yes) or 2D protofilaments (no).
- (Optional) Check the box to calculate terminal PF curvature. If checked, the terminal curvature (corresponding to the middle of the terminal dimer) will be reported in the command window of MATLAB and in term_angles.txt. It will also be recorded in the ang_vs_tip.txt file regardless of the box being checked.
- Specify the pixel size (nm).
- (Optional) Exclude protofilament tracings shorter than a given length from analysis by providing the lower length limit in tubulin dimers (default is 0 dimers).

Some optional steps before PF analysis:
- To manually browse PFs one by one: check the box to view individual PFs. Navigation: mouse-click to go forward, any keyboard key to go backward, "close all figures" to stop.
- To manually mark PFs in the 2D-mode: check the "mark individual PFs". Navigation: 1) first left mouse click on the PF to select, 2) second left mouse click to accept selection, right mouse click to deselect the protofilament, 3) press on scroll "close all figures" to stop.
- To exclude outlier protofilament tracings from analysis, create an 'excludePFs.txt' file in the folder containing your .dat files. The txt file should contain a vector of numbers for the protofilaments to be excluded, e.g., 1 20 21 30.

Press the "Start Processing" button.

Output txt files:
- all_angles - individual values of all angles between consecutive pairs of PF segments (deg/dimer).
- all_PF_lengths - individual values of PF lengths (nm).
- ang_vs_tip - dependence of average PF curvature on the distance from PF tip {distance(nm); curvature (deg/dimer); SD of curvature (deg/dimer)}.
- ang_vs_tip_fit - weighted linear fit of the dependence of average PF curvature on the distance from PF tip {distance(nm)}.
- Ang_ALL_mean_med_SD_N - average PF curvature (deg/dimer); median PF curvature (deg/dimer); SD of PF curvature (deg/dimer); N of angles.
- hist_ang_all - histogram of PF curvatures.
- hist_PF_length - histogram of PF lengths.
- MeanL_StdL_numL - mean PF length, SD pf PF length, N of PFs.
- term_angles - individual terminal curvature values.

Output jpeg/png files:
- angles - histogram of PF curvatures.
- histlength - histogram of PF lengths.
- plotPFs - all PF traces, smoothed.

--------------------------------------------------------------------
generate_test_PFs.m program
(This is an additional program, helpful in testing and troubleshooting the scripts above).
The program generates simulated PF tracings, having user-defined properties.

In the MATLAB Command Window:
- Run generate_test_PFs
- Enter the number of protofilaments.
- Enter the curvature, deg/dimer: 2.
- Enter the standard deviation of curvature, deg/dimer: 10.
- Enter the gradient of curvature, deg/dimer^2: 0.
- Select the type of PF length distribution: fixed / uniform / exponential / normal / gamma.
- Specify the parameters of the PF length distribution.
The program will create a folder named "simulated_coordinates" with "x.txt," "y.txt," and "z.txt" files containing PF coordinates.

