#==========================================
# +++ SIGNIFICANCE ASSESSMENT +++
# P-value estimation in a case/cont test
# by using mixture density with EM algorithm
# Exploration of the Likelihood space
# Gaussian or Beta mixture model
#==========================================


1) Requirements :

_ Distribution of loadings (standard output of Tax3 software)
_ loadings.csv file containing all loading observed before resampling procedure (standard output of Tax3 software)
_ List of variables

2) What does the analysis do ?

_ Estimation of the p-value in a case/cont 2 sided-test by using estimation of mixture density
_ Gaussian or Beta mixture density, respectively launched with Mixmod and R
_ Exploration of the entire space of Likelihood to take into consideration all local maximums of Likelihood
_ Mixture density graph is created into the results file
_ The runner_gaussian_mixure_debug_1.sh provides a confident interval (with SEM algorithm) for the isolated p-value
_ -> in this case, mixmod must be compiled with DEBUG=1 in the Mixmod directory (see <mixmod_dir>/LIB/MIXMOD/Util.h)

3) How to launch p-value estimation ?

_ 6 folders are required :
	* foler for datasets, loadings.csv and list of variable
	* folder for options files
	* folder for PGXIS R package (containing standard_plot_functions for graph)
	* folder for analysis (1 folder/variable)
	* folder for results (1 folder/variable)
	* folder for bash programs (go.sh, EM_functions.R, runner.sh with gaussian mixture or beta mixture)

_ launch the bash command "./<folder_for_bash_programs>/go.sh [-g or -b] options/options.dat"
_ see usage() bash function in the go.sh code for arguments
_ a log is created for each variable : analysis/<variable>/log_all.txt to follow the analysis



