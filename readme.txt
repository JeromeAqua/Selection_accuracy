1/ Need to create fake animal presence = synthetic tracks -- all in "Create_presence"
	a/ Download netcdf data for the environmental parameter organisms will be cueing on. Use Erdapp/ thredds /...
	b/ Run Biased_RW_* to create the biased random walks
	c/ When all the wanted environmental data have been run, run concatenate_RW.R to have everything in one neat "Fakeanimals.RData"

2/ Create animal pseudo-absences
	4 Different kinds of null models, all saved in all_tracks.RData
	(if want to plot, check Plot_absences folder)

3/ Extract environmental variables (satellite data + LCS) 
	a/ Extract_lcs_syntheticB.R will extract FTLE, Extract_sst_synthetic.R will extract both chl a and sst
	b/ Merge_lcs, Merge_chl, Merge_sst to merge all the datasets created in one. Data saved in Synthetic_all_LCS/sst/chl.RData. + Merge_all_with_0 to have the good datasets after having run the 0 case separately
	c/ Plot_histograms.R plots the density distribution of the resulting data

4/KS_tests
	a/ KS_chl, KS_LCS, and KS_sst commpute the KS tests and bootstrap confidence intervals. Note: We created new files (Synthetic_xxx_onlyxxx.RData) for all 3 variables to decrease the size of the dataset to open and save precious memory. Codes are run on clusters (biohen + caviness) to save time. Each file is usually launched for 1 value of kappa to have it easily parallelized. 
	b/ Bootstrap_plots plot the bootstrap estimates for the 3 variables




---------------------------------------------------------------------------
Re run for concentration = 0 to have true naive model at c = 0
1/ Run Biased_RW_SST0.R. Only changed l. 134, 135, 178, 230,289 compared to the normal files
!Use this one for SST, chl, and ftle. As there is no environmental selection they will be the same anyways.
2/ Create null models just for the concentration =0 files (concatenated after in all_tracks0.RData)
3/ Extract satellite data just for when kappa=0
