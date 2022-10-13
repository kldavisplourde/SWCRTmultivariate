The ZIP folder contains files for implementing the sample size and power calculations under a multivariate linear mixed model (MLMM) framework as shown in "Power analyses for stepped wedge designs with multivariate continuous outcomes" by Davis-Plourde, Taljaard, and Li.

List of Files:
1) Application.r = R file for generating the Application study including sensitivity analysis.
2) copri_test_HoopGir.r = R file for conducting the Simulation study.
3) EM_uncorrected_HoopGir.r = R file for EM Algorithm (used to fit MLMM).
4) gendata_copri_varCluster_HoopGir.r = R file for generating the Simulation study data.
5) Sim_Params.txt = Text file containing all variable inputs for each simulation scenario.
6) powerSampleCal_IU_HoopGir.r = R file for generating the power and sample size for the MLMM under a SW-CRT design using the intersection-union test.
7) powerSampleCal_omnibus_HoopGir.r = R file for generating the power and sample size for the MLMM under a SW-CRT design using the omnibus test..

NOTES:  1) This program requires the MASS package (comes preloaded into R so does not require installation). Simulation study additionally requires the installation of the lme4, doMC, doRNG, lmeInfo, mvtnorm, and numDeriv packages. Finally, the power and sample size function requires the installation of the mvtnorm package.
	2) You will need to change path names before running the programs. 
	3) Latest version of all files are available on GitHub: https://github.com/kldavisplourde/SWCRTmultivariate