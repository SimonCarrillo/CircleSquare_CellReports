PyCharm Analysis pipeline (Discard cells, AV/PTI, Isomap dimensionality reduction)

DiscardCells.py -> Discard cells by a given ratio and a given methodology (all cells kept, anti-cofiring removed, random removed)

Input: Temp_List_1000msAV_1000msKT_CircleSquare_….file (from cnmfe_analysis temporal.py)
Output: CellDiscarded_discarPop_Ratio_Temp_… .file  per every day per mouse, with information to serve as input for Isomap analysis.

main_Isomap.py -> Apply dimensionality reduction of the spiking data

Input: : CellDiscarded_discarPop_Ratio_Temp_… .file  per every day per mouse (output of DiscardCells.py)
Output: ISO_#dimensions_ CellDiscarded_ discarPop_ratio … . Mat file with a cell array. As many rows as number of days of recording (either 9 or 7). First column is raw_data input to IsoMap algorithm, second column is IsoMap reduced data, fifth column is the session identification (HMC, RCT, CYL).  Note that some codes need to number of dimensions removed from the filename to work!

tau_distance.py -> Computes the distance between pair of cells and its tau correlations 

Input: files from extraction pipeline (tracker videos, timestamps, excel spreadsheet, etc)
Output: Dist_M19.mat cell array with tau correlations and p-values for each pair of cells (first row), and distances between each pair of cells (second row).

Anticofiring_fields_analysis.py -> Calculates the number of anti-cofiring cell IDs overlapping between environments

Input: Temp_List_1000msAV_1000msKT_CircleSquare_….file (from cnmfe_analysis temporal.py)
Output: . Mat array with as many columns as days and 4 rows (row 1: overlap Env 1; row 2: overlap Env 2; row 3: overlap HMC; row 4: overlap Env 1 Env 2)


—————————————————————————————————————————————————
MATLAB Analysis pipeline (Manifolds, multi-stable angular, Kuiper’s statistic, best fit plane, Participation ratios)

manifolds2d_overlap.m -> Compute the proportion of overlap between 2d manifolds of CYL and RCT for a given mouse as well as the angle between best fit planes for the whole recording at once, 300 seconds. (n=6)

Input: ISO_CellDiscarded_kcRatio20TempCorr1000CircleSquare_M19_. Mat (output of main_Isomap.py)
Output: angle_M19_ratio8.mat and overlap_M19_ratio8.mat. As many rows as days in that mouse and three columns (col 1: all cells kept, col 2: anti-cofiring cells removed, col 3: random cells removed)

summary_overlap2dmanifolds.m -> Computes summary data in terms of mean and SEM for the overlap proportion between environment for the different ensembles (kc, neg, rand)  (n=6)

Input:angle_M19_ratio8.mat and overlap_M19_ratio8.mat for each animal (output of manifolds2d_overlap.m)
Output: Summary of statistics values that can be taken to Prism for final plotting

tau_dist.m -> fit a power law to analyze the scale-free properties of the tau correlations between a pair cell and its distance (n=6)

Input: Dist_M19.mat  (output of tau_distance.py)
Output: Power law fit of the distances between cells and their tau correlations, stats can be taken to Prism

bestfitplane_parametric.m -> Parametric analysis of the window size for best fit plane as a function of its residual (n=6) 

Input:  ISO_CellDiscarded_kcRatio20TempCorr1000CircleSquare_M19_. Mat (output of main_Isomap.py)
Output: Graph with window size for best fit plane vs residual for all the mice and its derivative

deltaangle_perrecording.m -> Computes angles between consecutive time steps for each recording per day of a single given mouse (n=6) 

Input:  ISO_CellDiscarded_kcRatio20TempCorr1000CircleSquare_M19_. Mat (output of main_Isomap.py)
Output: a list of .mat files of the following format, angles_all50.mat shifts_all50.mat, where the number corresponds to the angle threshold for that case.

angle_across_environments_kuipers.m ->  Computes the Kuiper test statistics for all the mice (n=6) from a range of angles between same an Dif environments

Input: ISO_CellDiscarded_kcRatio20TempCorr1000CircleSquare_M19_. Mat (output of main_Isomap.py)
Output: list of kuiper test statistics to be imported to Prism for analysis. Save in .mat file the angle between consecutive time steps for all the environments.

smoothness_localnessv2.m -> Compute smoothness and localness (n=6)

Input: ISO_CellDiscarded_kcRatio20TempCorr1000CircleSquare_M19_. Mat (output of main_Isomap.py)
Output: smoothness and localness values for all mice. Values can be taken to Prism for statistical analysis. Plots of the smoothness and localness with SEM shaded region.

partipation_ratio_analysis.m -> Compute participation ratios (using whole recordings, 300 s) for RCT, CYL, HMC as well as control groups for a single plane, plane shifting, Swiss roll and chance (n=6)

Input: ISO_CellDiscarded_kcRatio20TempCorr1000CircleSquare_M19_. Mat (output of main_Isomap.py)
Output: Participation ratio for every environment as matlab variables to be analyzed in Prism.

multistable_participation_ratio_isomapdim_parametric; multistable_participation_ratio_isomapdim_300s.m -> Compute Participation ratio for different Isomap dimension for all the mice (n=6)

Input: ISO_#dim_CellDiscarded_kcRatio2TempCorr1000CircleSquare_M29_.mat (output of main_isomap.py, by changing the number of dimensions)
Output: .mat file with participation ratio per environment, CYL, RCT, HMC with as many columns as number of Isomap dimensions investigated.

anti_cofiringpower.m -> Computes the anticofiring power for all the cells

Input: ISO_#dim_CellDiscarded_kcRatio2TempCorr1000CircleSquare_M29_.mat (output of main_isomap.py, by changing the number of dimensions)
Output: Histogram with distribution information of the anticofiring power across days and mice

———————————————————————————————————————
Functions needed for the main Matlab codes

Cellflat.m -> Flatten the tree dimensions of a given cell array to one row of cells

pwerfit.m -> Fit data using a given power law (used for analysis of scale-free tau correlations with distance)

circ_kuipertest.m; circ_samplecdf.m; kuipertable.mat -> Kuiper’s test statistical tool, mod from original to get K=nan if k<<

fitNormal.m -> Compute best fit plane and output the normal vector that defines the best fit plane

generate_data.m -> Generates a random set of points for defined geometries such as a Swiss roll, helix, etc

Stdshade.m -> Plotting function to create a shaded area (SEM) around the mean center line.

boxplotGroup.m -> To plot multiples box plots for grouped data (mainly two way ANOVAs), not used, all statistical analysis are conducted in Prism GraphPad.
