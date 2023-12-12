#!/bin/bash

# TODO :
# - allow for small variances : multiply input, divide output ...
#

#########################################################################################
#
#	EM mixture model
#        : selection of the best mixture of gaussian or beta distributions
#			: calculation of the pvalue
#
#
#	INPUTS
#         : observed loading
#         : loadings after case/control resampling
#
#  STATEGY
#     For each mixture model run:
#        - pre-run with a small number of loadings (500)
#        - Use obtained parameters for full run (usually 5000 loadings)
#        - SMALL_EM option
#        - Free models : Gaussian_pk_Lk_I , Beta_free_proportions
#
#     Main algorithm has 2 steps:
#        - step1 : starting from 1 distribution
#        - step2 : exploring around best model obtained from step1
#     Thorough Step1 only best for SNP level (small number of distributions expected)
#     Fast Step1 and thorough step2 best for gene level (high number of distributions expected)
#
#
#########################################################################################



verbose=false
vecho() {
if [[ $verbose == true ]] ; then
   echo $*
fi
}

# get minimum of chains to run, given number of distributions (used only if mode='adaptive')
getMinNChain() {
	ndist=$1

	if [ $ndist -lt 6 ] ; then
		n=`echo " ( $ndist * 3 ) " | bc`
	else
		n=`echo " ( $ndist * 2 ) " | bc`
	fi

	echo $n
}
# get minimum of chains to run, given number of distributions (used only if mode='adaptive')
getMaxNChain() {
	ndist=$1

	if [ $ndist -lt 6 ] ; then
		n=`echo " ( $ndist * 4 ) " | bc`
	else
		n=`echo " ( $ndist * 3 ) " | bc`
	fi

	echo $n
}

#---------------
# get parameters
#---------------

model=<MODEL>

# workdir
finaljobsdir=<JOBSDIR>
jobsdir=/ramdisk${finaljobsdir}

#
# results
#
finalhere=<HERE>
here=/ramdisk${finalhere}
mkdir -p $here
cp -R ${finalhere}/* ${here}

# where to find various functions
batch_dir=<BATCH_DIR>

# graphics
R_dir=<RDIR>

# variable being investigated
variable=<VAR>

# loading_type : component or projection
loading_type=<LOADING_TYPE>

# number of chains: fixed, adaptive
step1_chain_mode=<STEP1_CHAIN_MODE>
step2_chain_mode=<STEP2_CHAIN_MODE>

# Do we do a prerun with a small sample of the dataset to generate starting values?
# yes, no, alt (default)
# alt will do preruns 50% of the time as prerun with small sample set may bias results
doPreRun=<DO_PRE_RUN>

# PRERUN_INIT : SMALL_EM or RANDOM (default)
PreRunInitType=<PRE_RUN_INIT_TYPE>

# RunInitType ( if no PreRun or PreRun fails) : SMALL_EM (default) / RANDOM
RunInitType=<RUN_INIT_TYPE>

# STEP1
step1_n_chains=<STEP1_N_CHAINS>
step1_exploration_raised_limit=<STEP1_EXPLORATION_RAISED_LIMIT>
step1_exploration_higher_limit=<STEP1_EXPLORATION_HIGHER_LIMIT>
step1_max_dist=<STEP1_MAX_DIST>

# STEP 2
enable_step2=<ENABLE_STEP2>
step2_n_chains=<STEP2_N_CHAINS>
step2_exploration_raised_limit=<STEP2_EXPLORATION_RAISED_LIMIT>
step2_exploration_higher_limit=<STEP2_EXPLORATION_HIGHER_LIMIT>

# graph
dograph=<DOGRAPH>
title="<TITLE>"

PCA_component=<PCA_COMPONENT>

# number of crashes allowed
max_EM_error=<EM_ERROR>

# EM precision
EM_epsilon=<EM_EPSILON>
EM_nbiteration=<EM_NBITERATION>

# misc
minimize_IO=<MINIMIZE_IO>

# should we use all cpus to run the analysis? (yes / no) : not implemented yet
RunChainsInParallel=<RUN_CHAINS_IN_PARALLEL>

# Get pvalues CI
doCI=yes

### create folders of analysis ###
rm -rf $finaljobsdir"inputs" $finaljobsdir"jobs" $finaljobsdir"model_log" $finaljobsdir"log_files"
rm -rf $jobsdir"inputs" $jobsdir"jobs" $jobsdir"model_log" $jobsdir"log_files"
rm -rf $here"bic" $here"parameters"
mkdir -p $jobsdir"inputs" $jobsdir"jobs/R" $jobsdir"model_log" $jobsdir"log_files" $here"bic" $here"parameters/SEM" $here"parameters/step_1"

# quantile associated (loading observed)
quantile=`cut --delimiter=, --fields=2 ${here}loading.csv | tail -1`

# extraction of loadings after resamplings
cut --delimiter=, --fields=2 ${here}$variable.csv | tail -n+2 > ${jobsdir}inputs/mixture.dat

########################################################################################################
#
# DATA TRANFORM : MIRROR and FOLD
#
# As loadings have a non definite sign, it makes sense to make the set of distributed loading symetric around 0
# This is achived by duplicating the set with oposite values
# This is done practically line by line (so that it is still symetric locally : usefull for pre-run)
#
#
transform=<TRANSFORM>					# none, mirror or fold

mirror_symetry=<MIRROR_SYMETRY>		# make init parameters symetric
mirror_ndist_init=<MIRROR_NDIST_INIT>		# number of distributions to start from
mirror_ndist_step=<MIRROR_NDIST_STEP>		# step when increasing number of distributions
mirror_bic_correction=<MIRROR_BIC_CORRECTION>	# BIC correction (half the number of independent parameters)
########################################################################################################

# append to a file all its opposite values
# file has only one column and no header
# $1 : file name
mirror() {

	(
	cat <<- EOF
		x <- as.matrix(read.table(file="$1"))

		y <- array(dim=length(x)*2)

		for (i in 1:length(x) ) {
			y[i*2-1] <- x[i]
			y[i*2]   <- - x[i]
		}

		write(y, ncolumns=1, file="$1")
	EOF
	) > ${jobsdir}jobs/R/mirror.R

	( R --no-restore --no-save --no-readline -q < ${jobsdir}jobs/R/mirror.R ) &> ${jobsdir}jobs/R/log_mirror.txt

}

# make all values positive
# file has only one column and no header
# $1 : file name
fold() {

	(
	cat <<- EOF
		x <- as.matrix(read.table(file="$1"))

		y <- array(dim=length(x))

		for (i in 1:length(x) ) {
			y[i]   <- abs(x[i])

			# TODO : do it for both beta & gaussian for BIC consistency
			if ( ("$model" == "beta") && ( y[i]==0 ) ) {
				y[i] <- 0.00001 # or take half of lowest loading
			}

		}

		write(y, ncolumns=1, file="$1")
	EOF
	) > ${jobsdir}jobs/R/fold.R

	( R --no-restore --no-save --no-readline -q < ${jobsdir}jobs/R/fold.R ) &> ${jobsdir}jobs/R/log_fold.txt

}

if [[ $transform == "mirror" ]] ; then
   mirror ${jobsdir}inputs/mixture.dat
elif [[ $transform == "fold" ]] ; then
   fold ${jobsdir}inputs/mixture.dat
fi

########################################################################################################
#
# INTERNAL PARAMETERS
#
########################################################################################################

# get number of resamplings (including mirrored)
nb_line=`wc -l ${jobsdir}inputs/mixture.dat | cut --fields=1 --delimiter=' '`

# size of prerun : TODO : make this a user parameter
globalPrerunSize=$(echo "scale=0; $nb_line / 10 " | bc)

# maximum number of distributions
max_nbr_distributions=$(echo "scale=0;sqrt($nb_line)" | bc)


########################################################################################################
#
# GAUSSIAN AND BETA MIXTURE MODEL FUNCTIONS
#
########################################################################################################

# parameters
# -s : NbLines (size of dataset)
# -n : ListNbCluster (number of distribution)
# -e : StopRuleValue ( Epsilon )
# -i : StopRuleValue ( NbITerations )
# -d : DataFile
# -j : JobFile (.xem file)
# -l : LogFile
# -t : InitType ( SMALL_EM , USER)
# -f : InitFile
# -o : OutputPath (Beta only)

run_EM_model() {

	local NbLines
	local ListNbCluster
	local Epsilon
	local NbIteration
	local DataFile
	local JobFile
	local LogFile
	local InitType
	local InitFile
	local OutputPath

	# get options
	while getopts “s:n:e:i:d:j:l:t:f:o:” OPTION
	do
		case $OPTION in
		   s)
		       NbLines=$OPTARG
		       ;;
		   n)
		       ListNbCluster=$OPTARG
		       ;;
		   e)
		       Epsilon=$OPTARG
		       ;;
		   i)
		       NbIteration=$OPTARG
		       ;;
		   d)
		       DataFile=$OPTARG
		       ;;
		   j)
		       JobFile=$OPTARG
		       ;;
		   l)
		       LogFile=$OPTARG
		       ;;
		   t)
		       InitType=$OPTARG
		       ;;
		   f)
		       InitFile=$OPTARG
		       ;;
		   o)
		       OutputPath=$OPTARG
		       ;;
		  esac
	done

	if [[ "$model" == "gaussian" ]] ; then
		# mixmod command lines
		(
		cat <<- EOF
		NbLines
			${NbLines}
		PbDimension
			1
		NbCriterion
			1
		ListCriterion
			BIC
		NbNbCluster
			1
		ListNbCluster
			${ListNbCluster}
		NbModel
			1
		ListModel
			Gaussian_pk_Lk_I
		NbStrategy
			1
		InitType
			${InitType}
		EOF

		if [[ "$InitType" == "USER" ]]
		then
			cat <<- EOF
			InitFile
				${InitFile}
			EOF
		fi

		cat <<- EOF
		NbAlgorithm
			1
		Algorithm
			EM
		StopRule
			NBITERATION_EPSILON
		StopRuleValue
			${NbIteration} ${Epsilon}
		DataFile
			${DataFile}
		EOF
		) > $JobFile

		# launch mixmod and catch seg faults
		(
			mixmod $JobFile
			if [ $? = $(( 128 + 11 )) ]
			then
				echo -e "Analysis failed (MIXMOD ERROR: Segmentation fault)" 1>&2
			fi
		) 2>&1 > $LogFile

		# check exit status
		if grep "MIXMOD ERROR" $LogFile >/dev/null; then
			echo "error"
		else
			echo "ok"
		fi

	elif [[ "$model" == "beta" ]] ; then

		if [[ "$InitType" == "USER" ]] ; then
			SMALL_EM="FALSE"
		fi
		if [[ "$InitType" == "SMALL_EM" ]] ; then
			SMALL_EM="TRUE"
			InitFile="none"
		fi

		# TODO : implement RANDOM init & pass 'mirror' option to R

		# Launch beta/R
		(
		cat <<- EOF
		# The following files contains functions about the EM algorithm.
		source("${finaljobsdir}EM_functions.R")

		# dataset
		dataset=as.matrix(read.table("${DataFile}", header=FALSE))

		EOF

		if [[ $transform != "fold" ]] ; then
			cat <<- EOF
			# we get a dataset between 0 and 1 : this will require to correct the BIC by + 2*n*log(2)
			dataset=0.5*dataset+0.5

			EOF
		fi

		cat <<- EOF
		# Application of the EM algorithm
		EM_algorithm(dataset, ${ListNbCluster}, $Epsilon, $NbIteration, model="Beta_free_proportions", output.path="${OutputPath}", init.file="$InitFile", SMALL_EM=${SMALL_EM})
		EOF
		) > $JobFile

		# R
		(R --no-restore --no-save --no-readline -q < $JobFile) &>  $LogFile

		# check exit status
		if grep "Error" $LogFile >/dev/null; then
			echo "error"
		else
			echo "ok"
		fi


	fi


} # end run_EM_model function



########################################################################################################
#
# CHAIN FUNCTIONS
#
########################################################################################################


#====================================================================
# Runs several chains of mixture model
#
# Parameters :
# - $1 is the number of distributions
# - $2 is the minimal number of mixmod chains
# - $3 adaptive of fixed mode
#
# This function returns the number of the best chain according to the BIC criterion
# If all chains crash, it returns 0
#

launch_mixture()
{

# local variable i represent the number of distributions
local i=$1

# how many  chains ?
local chain

# ----------------------------------------------------
# determine number of chains to run

# how many independent distributions ?
local ndist
if [[ $transform == "mirror" ]] ; then
	temp=$1
	ndist=`echo " ( $temp + 1 ) / 2 " | bc`
else
	ndist=$1
fi

if [[ "$3" == "fixed" ]]
then
   chain=$2
fi

if [[ "$3" == "adaptive" ]]
then
   if [ $1 -eq 1 ]
   then
      # if one distribution : run always one chain
      chain=1
   else
   	nMinChain=`getMinNChain $ndist`
   	if [ $2 -gt $nMinChain ] ; then
   		nMinChain=$2
   	fi
      if [ $ndist -lt $nMinChain ]
      then
         # minimum number of chains
         chain=$nMinChain
      else
         # number of chains equals to number of (independant) distributions
         chain=$ndist
      fi
   fi
fi

# FIXME : infinite loop if chain=1
if [[ $chain -lt 2 ]] ; then
  chain=2
fi

# ----------------------------------------------------
# Local variables

# local array BICs represents the BICs for each chain
local BICs=()

# local array all_chains represents number of chains which succeed
local all_chains=()

# local index for vector
local k=0

# other local variables
local ic=0
local index=0
local pre_run=""
local main_run_error_counter=0
local pre_run_error_counter=0
local pre_defined=true
local InitType=""
local go=true
local count=0
local min=0
local best_chain=0

# get total number of data points
local nlines=`wc -l ${jobsdir}inputs/mixture.dat | cut --fields=1 --delimiter=' '`

# size of prerun
local prerunSize

# size of prerun when initiated
local prerunSizeInit

# increase size of prerun if many distributions
# that prevents EM from failing due to prerun data too sparse
if [[ $1 -gt 5 ]] ; then
	prerunSizeInit=`echo " (3 * $globalPrerunSize) / 2  " | bc`
elif [[ $1 -gt 10 ]] ; then
	prerunSizeInit=`echo " (4 * $globalPrerunSize) / 2  " | bc`
elif [[ $1 -gt 15 ]] ; then
	prerunSizeInit=`echo " (5 * $globalPrerunSize) / 2  " | bc`
elif [[ $1 -gt 20 ]] ; then
	prerunSizeInit=`echo " (6 * $globalPrerunSize) / 2  " | bc`
else
	prerunSizeInit=$globalPrerunSize
fi

# make folders for chain : jobs, model_log, log_files, BIC, parameters
mkdir -p ${jobsdir}jobs/dist_$i
mkdir -p ${jobsdir}model_log/dist_$i
mkdir -p ${jobsdir}log_files/dist_$i
mkdir -p ${here}bic/dist_$i
mkdir -p ${here}parameters/dist_$i
mkdir -p ${here}vars/dist_$i

#
# launch one chain nested function
#   take one parameter : the chain number
launch_one_chain() {

	local c=$1

	# log
   (
	echo -e "=============================================="
	echo -e "  Number of distributions : $i"
	echo -e "  Chain # $c "
	echo -e "=============================================="
   )  > ${jobsdir}log_files/dist_$i/log_chain_$c.txt

	# counter specific for main run
	main_run_error_counter=0

	go=true
	while [ $go == true ]
	do

		# log
		echo -e "\n########### Atempt $(( 1 + $main_run_error_counter)) ################" >> ${jobsdir}log_files/dist_$i/log_chain_$c.txt

		#--------------------------------
		# PRERUN
		#--------------------------------

		if [[ "${doPreRun}" == "yes" ]] ; then

			echo -e "\n------ PRE RUN -------" >> ${jobsdir}log_files/dist_$i/log_chain_$c.txt

			# counter specific for prerun
			pre_run_error_counter=0

			# variable which indicates whether the prerun is done or not
			pre_run="start"
			while [ "$pre_run" != "pre-run-done" ]
			do

				# ---------------------------------------------------------
				# Extraction of dataset with subset resamplings (for pre-run)

				# increase prerun size if main run or pre-run has alredy failed
				if [ $main_run_error_counter -eq 0 ] && [ $pre_run_error_counter -eq 0 ] ; then
					prerunSize=$prerunSizeInit
				else
					prerunSize=`echo " $prerunSizeInit + 100 * $main_run_error_counter + 100 * $pre_run_error_counter " | bc`
				fi

				# max limit for prerun
				if [ $prerunSize -gt $nlines ] ; then
					prerunSize=$nlines
				fi

				# extract data sample
				(
				cat <<- EOF
					x <- as.matrix(read.table(file="${jobsdir}inputs/mixture.dat"))

					y <- sample(x, size=$prerunSize )

					write(y, ncolumns=1, file="${jobsdir}inputs/data_$c.dat")
				EOF
				) > ${jobsdir}jobs/R/sample.R

				( R --no-restore --no-save --no-readline -q < ${jobsdir}jobs/R/sample.R ) &> ${jobsdir}jobs/R/log_sample.txt

				# log
				(
				echo -e "\nPre-run :"
				echo -e "   - prerunsize = ${prerunSize}"
				echo -e "   - nlines     = ${nlines}"
				) >> ${jobsdir}log_files/dist_$i/log_chain_$c.txt

				#-----------------------------------------------------------------------------
				# PRERUN EM
				# TODO : tell R to start with mirrored values if necessary

				echo -e "\nPre-run in progress ..." >> ${jobsdir}log_files/dist_$i/log_chain_$c.txt

				# call model
				# -s : NbLines (size of dataset)
				# -n : ListNbCluster (number of distribution)
				# -e : StopRuleValue ( Epsilon )
				# -i : StopRuleValue ( NbIteration )
				# -d : DataFile
				# -j : JobFile (.xem file)
				# -l : LogFile
				# -t : InitType ( SMALL_EM , USER)
				# -f : InitFile
				# -o : Output File Path (Beta only)

				# beware : !!! no trailing tabs after backslashes !!!
				if [[ "$model" == "gaussian" ]] ; then
					exit_status=$(run_EM_model 	\
									-s ${prerunSize} 																		\
									-n $i																						\
									-e ${EM_epsilon}																		\
									-i ${EM_nbiteration}																	\
									-d "${jobsdir}inputs/data_$c.dat"												\
									-j "${jobsdir}jobs/dist_$i/mixmod_options_prerun_chain_$c.xem"			\
									-l "${jobsdir}model_log/dist_$i/_log_prerun_chain_$c.txt"				\
									-t ${PreRunInitType} )

				elif [[ "$model" == "beta" ]] ; then
					exit_status=$(run_EM_model 	\
									-s ${prerunSize} 																		\
									-n $i																						\
									-e ${EM_epsilon}																		\
									-i ${EM_nbiteration}																	\
									-d "${jobsdir}inputs/data_$c.dat"												\
									-j "${jobsdir}jobs/dist_$i/EM_options_prerun_chain_$c.R"					\
									-l "${jobsdir}model_log/dist_$i/_log_prerun_chain_$c.txt"				\
									-o "${here}"																			\
									-t ${PreRunInitType} )
				fi

				# check exit status
				if [[ "$exit_status" != "ok" ]] ; then

					# Error
					(( pre_run_error_counter+=1 ))
					# log
					echo -e "\nPre-run $pre_run_error_counter failed" >> ${jobsdir}log_files/dist_$i/log_chain_$c.txt
					if [ $pre_run_error_counter -eq $max_EM_error ]
					then
						echo -e "\nMax pre-run limit reached" >> ${jobsdir}log_files/dist_$i/log_chain_$c.txt
						pre_run="pre-run-done"
						pre_defined=false

					fi

				else	# sucess
					# if there is no error, then we retreive model parameters

					# make init parameters symetric if required
					if [[ $transform == "mirror" ]] && [[ $mirror_symetry == true ]] ; then

						echo -e "\nMaking distributions symetric" >> ${jobsdir}log_files/dist_$i/log_chain_$c.txt

						if [[ "$model" == "gaussian" ]] ; then
							# TODO : make sure we always have a dist centered on zero if odd number of dist
							# keep distributions having a positive (or zero) mean and duplicate
							(
							cat <<- EOF
								# import parameters
								parameters = as.matrix(read.table(file="BICparameter.txt"))

								# how many distributions ?
								K = nrow(parameters)/3
								index=1:( K * 3)

								# proportions, means & variances
								p = as.vector(parameters[which(index%%3 == 1),1])
								m = as.vector(parameters[which(index%%3 == 2),1])
								v = as.vector(parameters[which(index%%3 == 0),1])

								# new parameters
								param2 = array(dim=length(parameters))

								# TODO : if odd K, force lowest mean to zero

								j = 1
								for (i in 1:K) {
									if ( (m[i] == 0) && (j<nrow(parameters)) ) {
										# keep as is if mean is zero
										param2[j]   = p[i]
										param2[j+1] = m[i]
										param2[j+2] = v[i]
										j = j + 3
									} else if ( (m[i] > 0) && (j<nrow(parameters)) ) {
										# keep if mean positive
										param2[j]   = p[i]
										param2[j+1] = m[i]
										param2[j+2] = v[i]
										# and duplicate with mirror
										param2[j+3] = p[i]
										param2[j+4] = - m[i]
										param2[j+5] = v[i]
										j = j + 6
									}
								}

								if (j < nrow(parameters)) {
									# not enough distributions - try the opposite direction
									j=1
									for (i in 1:K) {
										if ( (m[i] < 0 ) && (j<nrow(parameters)) ) {
											# keep if mean negative
											param2[j]   = p[i]
											param2[j+1] = m[i]
											param2[j+2] = v[i]
											# and duplicate with mirror
											param2[j+3] = p[i]
											param2[j+4] = - m[i]
											param2[j+5] = v[i]
											j = j + 6
										}
									}
								}

								# bring sum of proportions to 1
								s = sum(as.vector(param2[which(index%%3 == 1)]))
								param2[which(index%%3 == 1)] = param2[which(index%%3 == 1)] / s

								write(param2, ncolumns=1, file="${here}parameters/dist_$i/parameter_prerun_chain_$c.init")
								write(parameters, ncolumns=1, file="${here}parameters/dist_$i/parameter_prerun_premirror_chain_$c.init")
							EOF
							) > ${jobsdir}jobs/R/param_${i}_${c}.R

						elif [[ "$model" == "beta" ]] ; then

							# keep positive distributions and duplicate
							(
							cat <<- EOF
								# import parameters
								parameters = as.matrix(read.table(file="${here}parameters.txt"))

								# how many distributions ?
								K = nrow(parameters)/3
								index=1:( K * 3)

								# proportions, means & variances
								p = as.vector(parameters[which(index%%3 == 1),1])
								a = as.vector(parameters[which(index%%3 == 2),1])
								b = as.vector(parameters[which(index%%3 == 0),1])

								# new parameters
								param2 = array(dim=length(parameters))

								# TODO : if odd K, force closest a,b to their average

								j = 1
								for (i in 1:K) {
									if ( (a[i] == b[i]) && (j<nrow(parameters)) ) {
										# keep as is if mean is 0.5
										param2[j]   = p[i]
										param2[j+1] = a[i]
										param2[j+2] = b[i]
										j = j + 3
									} else if ( (a[i] > b[i]) && (j<nrow(parameters)) ) {
										# keep if mean = a/(a+b) > 0.5
										param2[j]   = p[i]
										param2[j+1] = a[i]
										param2[j+2] = b[i]
										# and duplicate with mirror
										param2[j+3] = p[i]
										param2[j+4] = b[i]
										param2[j+5] = a[i]
										j = j + 6
									}
								}

								if (j < nrow(parameters)) {
									# try other side
									j = 1
									for (i in 1:K) {
										if ( (b[i] > a[i]) && (j<nrow(parameters)) ) {
											# keep if mean = a/(a+b) < 0.5
											param2[j]   = p[i]
											param2[j+1] = a[i]
											param2[j+2] = b[i]
											# and duplicate with mirror
											param2[j+3] = p[i]
											param2[j+4] = b[i]
											param2[j+5] = a[i]
											j = j + 6
										}
									}
								}

								# bring sum of proportions to 1
								s = sum(as.vector(param2[which(index%%3 == 1)]))
								param2[which(index%%3 == 1)] = param2[which(index%%3 == 1)] / s

								write(param2, ncolumns=1, file="${here}parameters/dist_$i/parameter_prerun_chain_$c.init")
								write(parameters, ncolumns=1, file="${here}parameters/dist_$i/parameter_prerun_premirror_chain_$c.init")
							EOF
							) > ${jobsdir}jobs/R/param_${i}_${c}.R

						fi # model==beta

						( R --no-restore --no-save --no-readline -q < ${jobsdir}jobs/R/param_${i}_${c}.R ) &> ${jobsdir}jobs/R/log_param_${i}_${c}.txt

					else
						# symetric parameters not required : copy all parameters as is
						if [[ "$model" == "gaussian" ]] ; then
							cp BICparameter.txt ${here}parameters/dist_$i/parameter_prerun_chain_$c.init
						elif [[ "$model" == "beta" ]] ; then
							cp ${here}parameters.txt ${here}parameters/dist_$i/parameter_prerun_chain_$c.init
						fi

					fi # symetry

					# the pre-run with 500 resamplings is done and succeed
					pre_run="pre-run-done"
					echo -e "\nPre-run successfull" >> ${jobsdir}log_files/dist_$i/log_chain_$c.txt

					# we will start EM with pre-defined values
					pre_defined=true

				fi

			done	# end of pre-run

		fi # ${doPreRun} == true

		#--------------------------------
		# START EM WITH ALL RESAMPLINGS
		#--------------------------------

		echo -e "\n------ MAIN ANALYSIS -------" >> ${jobsdir}log_files/dist_$i/log_chain_$c.txt

		# check if the pre-run succeed or not
		if [[ $pre_defined == true ]] && [[ "${doPreRun}" == "yes" ]] ; then
			echo -e "\nInitialization with pre-run output model parameters" >> ${jobsdir}log_files/dist_$i/log_chain_$c.txt
			# initiation with pre-defined values
			InitType=USER
		else
			if [[ $transform != "mirror" ]] || [[ $mirror_symetry != true ]] ; then
				# SMALL_EM option for initialization
				echo -e "\nInitialization with $RunInitType strategy" >> ${jobsdir}log_files/dist_$i/log_chain_$c.txt
				InitType=$RunInitType
			else
				# symetric random starting values if mirror
				echo -e "\nInitialization with random +mirrored+ values" >> ${jobsdir}log_files/dist_$i/log_chain_$c.txt
				InitType=USER
				if [[ "$model" == "gaussian" ]] ; then
					# generate
					(
					cat <<- EOF
						# number of distributions to mirror
						K = ceiling($i / 2) ;

						set.seed(`date +"%s"`)

						# parameters
						param2 = array(dim=3*$i)

						j=1
						if ( $i > 1) {
							for (i in 1:floor( $i / 2) ) {
								p = runif(1, min=0.05, max=1)
								m = runif(1, min=0.02, max=0.5)
								v = runif(1, min=0.001, max=0.01)
								param2[j]   = p
								param2[j+1] = m
								param2[j+2] = v
								param2[j+3] = p
								param2[j+4] = - m
								param2[j+5] = v
								j = j + 6
							}
						}

						# if odd number of dist, add a dist centered on zero
						if (K != floor( $i / 2 ) ) {
							p = runif(1, min=0.05, max=1)
							m = 0
							v = runif(1, min=0.001, max=0.01)
							param2[j]   = p
							param2[j+1] = m
							param2[j+2] = v
						}

						index=1:( $i * 3)

						s = sum(as.vector(param2[which(index%%3 == 1)]))
						param2[which(index%%3 == 1)] = param2[which(index%%3 == 1)] / s

						write(param2, ncolumns=1, file="${here}parameters/dist_$i/parameter_prerun_chain_$c.init")
					EOF
					) > ${jobsdir}jobs/R/mirror_param.R

				elif [[ "$model" == "beta" ]] ; then

					# generate
					(
					cat <<- EOF
						# number of distributions to mirror
						K = ceiling($i / 2) ;

						set.seed(`date +"%s"`)

						# parameters
						param2 = array(dim=3*$i)

						j=1
						if ( $i > 1) {
							for (i in 1:floor( $i / 2) ) {
								p = runif(1, min=0.05, max=1)
								a = runif(1, min=0.1, max=10000)
								b = runif(1, min=0.1, max=10000)
								param2[j]   = p
								param2[j+1] = a
								param2[j+2] = b
								param2[j+3] = p
								param2[j+4] = b
								param2[j+5] = a
								j = j + 6
							}
						}

						# if odd number of dist, add a dist centered on zero
						if (K != floor( $i / 2 ) ) {
							p = runif(1, min=0.05, max=1)
							a = runif(1, min=0.1, max=10000)
							param2[j]   = p
							param2[j+1] = a
							param2[j+2] = a
						}

						index=1:( $i * 3)

						s = sum(as.vector(param2[which(index%%3 == 1)]))
						param2[which(index%%3 == 1)] = param2[which(index%%3 == 1)] / s

						write(param2, ncolumns=1, file="${here}parameters/dist_$i/parameter_prerun_chain_$c.init")
					EOF
					) > ${jobsdir}jobs/R/mirror_param.R

				fi

			   ( R --no-restore --no-save --no-readline -q < ${jobsdir}jobs/R/mirror_param.R ) &> ${jobsdir}jobs/R/log_mirror_param_R.txt

		   fi # mirror

		fi # pre_defined


		# call model
		# -s : NbLines (size of dataset)
		# -n : ListNbCluster (number of distribution)
		# -e : StopRuleValue ( Epsilon )
		# -i : StopRuleValue ( NbIteration )
		# -d : DataFile
		# -j : JobFile (.xem file)
		# -l : LogFile
		# -t : InitType ( SMALL_EM , USER)
		# -f : InitFile
		# -o : Output File Path (Beta only)

		# beware : !!! no trailing tabs after backslashes !!!
		if [[ "$model" == "gaussian" ]] ; then
			exit_status=$(run_EM_model 	\
							-s ${nb_line} 																			\
							-n $i																						\
							-e ${EM_epsilon}																		\
							-i ${EM_nbiteration}																	\
							-d "${jobsdir}inputs/mixture.dat"												\
							-j "${jobsdir}jobs/dist_$i/mixmod_options_chain_$c.xem"					\
							-l "${jobsdir}model_log/dist_$i/_log_chain_$c.txt"							\
							-t ${InitType}																			\
							-f "${here}parameters/dist_$i/parameter_prerun_chain_$c.init"	)
		elif [[ "$model" == "beta" ]] ; then
			exit_status=$(run_EM_model 	\
							-s ${nb_line} 																			\
							-n $i																						\
							-e ${EM_epsilon}																		\
							-i ${EM_nbiteration}																	\
							-d "${jobsdir}inputs/mixture.dat"												\
							-j "${jobsdir}jobs/dist_$i/EM_options_chain_$c.R"							\
							-l "${jobsdir}model_log/dist_$i/_log_chain_$c.txt"							\
							-o "${here}"																			\
							-t ${InitType}																			\
							-f "${here}parameters/dist_$i/parameter_prerun_chain_$c.init"	)
		fi

		# check exit status
		if [[ "$exit_status" != "ok" ]] ; then
			# EM crashed
			(( main_run_error_counter+=1 ))
			echo -e "Analysis #$main_run_error_counter failed : re trying !!!\n\n" >> ${jobsdir}log_files/dist_$i/log_chain_$c.txt
			if [ $main_run_error_counter -eq $max_EM_error ] ; then
				# panic if too many errors
				echo -e "\n\nPanic : All analyses failed - Ending\n" >> ${jobsdir}log_files/dist_$i/log_chain_$c.txt
				go=false
			fi
		else
			# EM successfull
			echo -e "\nAnalysis successfull" >> ${jobsdir}log_files/dist_$i/log_chain_$c.txt

			# get results : we just need BIC and output model parameters
			if [[ "$model" == "gaussian" ]] ; then
				# Save of mixmod outputs into directory : BIC criterion, parameters
				cp BICparameter.txt ${here}parameters/dist_$i/parameter_chain_$c.txt
				cp BICnumericStandard.txt ${here}bic/dist_$i/bic_chain_$c.txt

				# get BIC - line 7 of 'BICnumericStandard.txt'
				bic=$(sed -n 7p BICnumericStandard.txt)

			elif [[ "$model" == "beta" ]] ; then

				# Save R outputs : BIC criterion, parameters
				cp ${here}parameters.txt ${here}parameters/dist_$i/parameter_chain_$c.txt
				cp ${here}BIC.txt        ${here}bic/dist_$i/bic_chain_$c.txt

				# get BIC
				bic=$(sed -n 1p ${here}BIC.txt)
			fi

			# BIC CORRECTION related to dataset*0.5+0.5
			if [[ "$model" == "beta" ]] && [[ $transform != "fold" ]] ; then
				bic=$(echo "scale=20 ; $bic + 2*$nb_line*l(2) " | bc -l )
			fi

			# Correct BIC if mirror : half the number of independent parameters
			# BIC= -2.L + nb_free_param*ln(n)
			# ici nb_free_param=(3*K-1)
			if [[ $transform == "mirror" ]] && [[ $mirror_bic_correction == true ]] ; then
				# This corrects for the number of independent distributions
				# Gives bad convergence (too many distributions)
   			# BICs[k]=`echo "scale=20 ; $bic - (3*$i-1)*l($nb_line) + (3*$i/2-1)*l($nb_line) " | bc -l`

   			# This corrects for :
   			#    i) non independance of dataset : LogL divided by 2
   			#   ii) even number of dependent distributions : 3K/2 -1
   			#   iii) odd number of dependent distributions (i.e. one distrib with a mean=0) : 3.(K-1)/2 +2 -1
   			LogL=$(echo "scale=20; -0.5 * ( $bic - (3*$i-1)*l($nb_line) ) " | bc -l )
   			if [[ $(( $i % 2 )) -eq 0 ]] ; then
   				# even
   				BICs[k]=$(echo "scale=20 ; -1 * ($LogL) + (3*$i/2-1)*l( $nb_line / 2) " | bc -l )
   			else
   				# odd
   				BICs[k]=$(echo "scale=20 ; -1 * ($LogL) + (3*($i-1)/2 + 2 - 1 )*l( $nb_line / 2) " | bc -l )
				fi
   		else
   		   BICs[k]=$bic
   		fi


   		# store bic
   		echo ${BICs[k]} > ${here}bic/dist_$i/bic_final_chain_$c.txt

			# lock to store all these results with no conflicts
			lockfile -1 ${here}vars/dist_$i/source_chain.txt.lock

			# get pvalue for these parameters
			if [[ "$model" == "gaussian" ]] ; then
				(
				cat <<- EOF
					# import parameters
					parameters=as.matrix(read.table("${here}parameters/dist_$i/parameter_chain_$c.txt"))

					# how many gaussian distributions ?
					K=nrow(parameters)/3 ; index=1:(3*K) ;

					# proportions, means & variances
					proportions=as.vector(parameters[which(index%%3 == 1),1])
					proportions=proportions  + (1-sum(proportions))/length(proportions)
					means=as.vector(parameters[which(index%%3 == 2),1])
					variances=as.vector(parameters[which(index%%3 == 0),1])

					# pvalue
					u = sum(proportions*pnorm(abs($quantile),means,sqrt(variances),lower.tail=FALSE)) # right
					v = sum(proportions*pnorm(-abs($quantile),means,sqrt(variances),lower.tail=TRUE)) # left

					w = sum(proportions*pnorm(1,means,sqrt(variances),lower.tail=FALSE)) # excess on the right
					z = sum(proportions*pnorm(-1,means,sqrt(variances),lower.tail=TRUE)) # excess on the left

					if ( "$transform" == "fold" ) {
						pvalue=u
						pvalue_excess_over1=w
					} else {
						pvalue=u+v
						pvalue_excess_over1=w+z
					}

					# save results
					results=matrix(0,ncol=2,nrow=1)
					results[1,1]=pvalue
					results[1,2]=pvalue_excess_over1

					colnames(results)=c("pvalue","pvalue_excess_over1")
					write.csv(results, file="${here}parameters/dist_$i/pvalues_chain_$c.csv", quote=FALSE, row.names=FALSE)
				EOF
				) > ${jobsdir}jobs/R/pvalue_$c.R

			elif [[ "$model" == "beta" ]] ; then

				(
				cat <<- EOF
				# import parameters
				parameters=as.matrix(read.table("${here}parameters/dist_$i/parameter_chain_$c.txt"))

				# how many distributions ?
				K=nrow(parameters)/3 ; index=1:(3*K) ;

				# proportions, means & variances
				proportions=as.vector(parameters[which(index%%3 == 1),1])
				proportions=proportions  + (1-sum(proportions))/length(proportions)
				means=as.vector(parameters[which(index%%3 == 2),1])
				variances=as.vector(parameters[which(index%%3 == 0),1])

				# proportions, alpha and beta
				proportions=as.vector(parameters[which(index%%3 == 1),1])
				proportions=proportions + (1-sum(proportions))/length(proportions)
				alpha=as.vector(parameters[which(index%%3 == 2),1])
				beta=as.vector(parameters[which(index%%3 == 0),1])

				# pvalue
				EOF
				if [[ $transform == "fold" ]] ; then
					cat <<- EOF
					u = sum(proportions*pbeta(abs($quantile),alpha,beta,lower.tail=FALSE)) 	# right
					pvalue=u
					EOF
				else
					cat <<- EOF
					u = sum(proportions*pbeta(0.5+0.5*abs($quantile),alpha,beta,lower.tail=FALSE)) 	# right
					v = sum(proportions*pbeta(0.5-0.5*abs($quantile),alpha,beta,lower.tail=TRUE))		# left
					pvalue=u+v
					EOF
				fi
				cat <<- EOF

				# save results
				results=matrix(0,ncol=2,nrow=1)
				results[1,1]=pvalue
				results[1,2]=0 # excess over one for compatibility with gaussian function

				colnames(results)=c("pvalue","pvalue_excess_over1")
				write.csv(results, file="${here}parameters/dist_$i/pvalues_chain_$c.csv", quote=FALSE, row.names=FALSE)

				EOF
				) > ${jobsdir}jobs/R/pvalue_$c.R

			fi

			( R --no-restore --no-save --no-readline -q < ${jobsdir}jobs/R/pvalue_$c.R ) &> ${jobsdir}jobs/R/log_pvalue_${c}_R.txt

         # LOG
			(
				if [ -f ${here}parameters/dist_$i/pvalues_chain_$c.csv ] ; then
			      echo -n "$variable,$PCA_component,$model,$i,$c,${BICs[k]},"
			      tail -1 ${here}parameters/dist_$i/pvalues_chain_$c.csv
				fi
			) >> ${here}"bics.csv"

			# which chain is it ?
			all_chains[k]=$c

			# store vars
			echo BICs['$k']=${BICs[k]} >> ${here}vars/dist_$i/source_chain.txt
			echo all_chains['$k']=${all_chains[k]} >> ${here}vars/dist_$i/source_chain.txt
			echo '(( k+=1 ))' >> ${here}vars/dist_$i/source_chain.txt

			# remove lock
			rm -f ${here}vars/dist_$i/source_chain.txt.lock

			# success : increment k
			(( k+=1 ))

			# stop this chain
			go=false

		fi

	# end of this chain
	done


} # end of launch_one_chain nested function

#-----------------
# Loop over chains
#-----------------

# not implemented properly
if [[ x"${RunChainsInParallel}" == x"yes_test" ]] ; then
	# get number of online cores
	nthreads=`cat /sys/devices/system/cpu/cpu*/online | paste -sd+ | bc`
	# add 1 for cpu0 (always online)
	(( nthreads+=1 ))
else
	nthreads=1
fi

# init
ic=1
runMoreChain=1
rm -f ${here}vars/dist_$i/source_chain.txt
rm -f ${here}vars/dist_$i/source_chain.txt.lock
echo 'k=0' > ${here}vars/dist_$i/source_chain.txt

# loop
while [  ${runMoreChain} -eq 1 ]
do
	while [  ${ic} -le ${chain} ]
	do

		for (( t=1; t<=nthreads; t++ ))
		do
			launch_one_chain ${ic} &
			(( ic+=1 ))
		done

		wait

	# end of loop overs chains
	done

	source ${here}vars/dist_$i/source_chain.txt

	# is there enough redundancy in BICs?
	# if not do one more chain
	if [[ $k -ge 2 ]] ; then
		rm ${jobsdir}log_files/dist_$i/bics_chain_${ic}.dat 2> /dev/null

		for (( count=0; count<$k; count++ ))
		do
			echo ${BICs[$count]} >> ${jobsdir}log_files/dist_$i/bics_chain_${ic}.dat
		done

		(
		cat <<- EOF

			#
			# we check that the lowest BICs are redundant
			#

			dat = as.matrix(read.table(file="${jobsdir}log_files/dist_$i/bics_chain_${ic}.dat"))

			# make it more or less stringent here :
			dat2 = round(dat)
			# dat2 = round(dat*10)

			# number of occurence of lowest BIC
			n <- length( dat2[dat2[]==min(dat2)] )

			write(n, ncolumns=1, file="${jobsdir}log_files/dist_$i/bics_nbest_chain_${ic}.dat")

		EOF
		) > ${jobsdir}jobs/R/unique.R

		( R --no-restore --no-save --no-readline -q < ${jobsdir}jobs/R/unique.R ) &> ${jobsdir}jobs/R/log_unique_R.txt

		# solo BICs : should be zero
		nbestBICs=`cat ${jobsdir}log_files/dist_$i/bics_nbest_chain_${ic}.dat`

		# test : lowest bic should be present twice
		echo -e "\n\nAll chains done - Testing for redundancy:" >> ${jobsdir}log_files/dist_$i/log_chain_${ic}.txt
		if [ $nbestBICs -lt 2 ] ; then
			if [ $chain -gt `getMaxNChain $i` ] ; then
				echo "  : Not enough BIC redundancy, but too many chains done" >> ${jobsdir}log_files/dist_$i/log_chain_${ic}.txt
				runMoreChain=0
			else
				echo "  : Not enough BIC redundancy : doing more chains" >> ${jobsdir}log_files/dist_$i/log_chain_${ic}.txt
				(( chain+=1 ))
				runMoreChain=1
			fi
		else
			echo "  : Enough BIC redundancy" >> ${jobsdir}log_files/dist_$i/log_chain_${ic}.txt
			runMoreChain=0
		fi

	fi


# end of RunMoreChain
done


# cat ${here}vars/dist_$i/source_chain.txt > /volatile/source_chain_dist$i.txt

#------------------------------------------------------------------------
# All chains are done : Selection of the best chain according to the BIC
#------------------------------------------------------------------------

echo "All chains are done : Selection of the best chain according to the BIC" >> ${jobsdir}log_files/dist_$i/log_chains.txt

# check if at least one chain succeed
if [ $k -ge 1 ]
then
	# check if at least two chains succeed
	if [ $k -ge 2 ]
	then
		# start values for the best chain
		count=1
		min=${BICs[0]}
		best_chain=${all_chains[0]}

		# Search the minimum
		while [ $count -lt $k ]
		do
			if [ "$(echo "$min > ${BICs[count]}" | bc)" -eq 1 ]
			then
				# min of BIC for this chain
				min=${BICs[count]}
				best_chain=${all_chains[count]}
			fi
			(( count+=1 ))
		done
	else
		min=${BICs[0]}
		best_chain=${all_chains[0]}
	fi
else
	best_chain=0
	min=0
fi

echo "   - Number of successful chains : $k " >> ${jobsdir}log_files/dist_$i/log_chains.txt
echo "   - Best chain : $best_chain " >> ${jobsdir}log_files/dist_$i/log_chains.txt
echo "   - min BIC    : $min " >> ${jobsdir}log_files/dist_$i/log_chains.txt

### END OF EM FUNCTION : return the number of the best chain
echo $best_chain

} # end of EM function


########################################################################################################
#
# MAIN ALGORITHM
#
########################################################################################################

# Run from the analysis directory
cd $jobsdir

(
# Calculation time : beginning
begin=$(date +%s)

echo -e "\n==============================================="
echo -e "\n 	PGXIS - SIGNIFICANCE ASSESSMENT"
echo -e "\n 	EM algorithm with $model mixture model"
echo -e "\n     - variable      : $variable"
echo -e "     - PCA_component : $PCA_component"
echo -e "     - title         : "$title
echo -e "\n==============================================="


#===============================================================
#  STEP 1 : initial exploration
#  - allows thorough exploration : then step2 not needed
#  - allows fast exploration : then step2 needed
#===============================================================

echo -e "\n----------------------------------------------------------"
echo -e "STEP 1"
echo -e "----------------------------------------------------------"

# log parameters
echo -ne "\nParameters"
echo -ne "\n  data tranform            : $transform"
if [[ $transform == "mirror" ]] ; then
	echo -ne "\n     Mirror parameters:"
	echo -ne "\n       symetric init    : $mirror_symetry"
	echo -ne "\n       initial nb dist  : $mirror_ndist_init"
	echo -ne "\n       step for nb dist : $mirror_ndist_step"
	echo -ne "\n       BIC correction   : $mirror_bic_correction"
fi
echo -ne "\n  number of resamplings    : $nb_line"
if [[ $transform == "mirror" ]] ; then
   echo -ne " (including mirrored resamplings)"
fi
echo -ne "\n  model                    : $model"
echo -ne "\n  number of chains         : $step1_chain_mode , minimal number=$step1_n_chains"
echo -ne "\n  exploration raised limit : $step1_exploration_raised_limit"
echo -ne "\n  exploration higher limit : $step1_exploration_higher_limit"
echo -ne "\n  exploration ndist  limit : $step1_max_dist"
echo -ne "\n  max number of crashes    : $max_EM_error"
echo -ne "\n  EM epsilon               : $EM_epsilon"
echo -ne "\n  EM nb iteration          : $EM_nbiteration"
echo -ne "\n  doPreRun                 : $doPreRun"
echo -ne "\n  PreRunInitType           : $PreRunInitType"
echo -ne "\n  RunInitType              : $RunInitType"
echo -ne "\n\n"

# i represents the best number of distributions, i.e. providing the lowest BIC : BIC_star

# starting point and progression step according to mirror
if [[ $transform == "mirror" ]] ; then
   # starting with two distribution
   ndist_init=$mirror_ndist_init
   # even number of distributions
   ndist_step=$mirror_ndist_step
else
   # starting with one distribution
   ndist_init=1
   # any distributions
   ndist_step=1
fi
ndist=$ndist_init

# number of times bic higher
nhigher=0

# number of times bic has raised
nraised=0

# used if data mirrored and BIC corrected
nraised_odd=0   # for an odd number of distribution
nraised_even=0  # for an even number of distribution

# initialize BICs
previous_BIC_tilde=10^100
previous_evenBIC_tilde=10^100
previous_oddBIC_tilde=10^100

# any success
anysuccess=false


# log headers
echo "ndist,bic" > ${here}bic/recap_step_1.csv
if [[ $loading_type == "component" ]] ; then
	echo "variable,PCA_component,model,n_distributions,id_chain,bic,pvalue,pvalue_excess_over1" > ${here}"bics.csv"
else
	echo "variable,PCA_dimension,model,n_distributions,id_chain,bic,pvalue,pvalue_excess_over1" > ${here}"bics.csv"
fi


# loop over distributions
while [ $ndist -le $max_nbr_distributions ]
do

   # Allows several tries with more chains and more gaussians if too many crashes
   success=0
   tries=0	# number of chain crashes
   nraise=0
   current_ndist=$ndist
   ncrash=0 # number of successive distribution crashes
   while [ $success -eq 0 ] && [ $current_ndist -le $max_nbr_distributions ] && [ $ncrash -le $max_EM_error ]
   do
      # log
      if [ $current_ndist -eq 1 ]
      then
         echo -en "Running $current_ndist distribution  ... "
      else
         echo -en "Running $current_ndist distributions ... "
      fi

      # launch mixture model
      local_n_chains=`echo "$step1_n_chains + $tries" | bc`
   	c=`launch_mixture $current_ndist $local_n_chains $step1_chain_mode`
   	(( tries+=1 ))

   	# test outcome
   	if [[ "$c" != "0" ]]
   	then
   	   # success
         success=1
         ndist=$current_ndist
      else
         # failure : rerun with one more chain
         #           unless too many errors
         #             : in this case increase number of distributions
         if [ $tries -ge $max_EM_error ]
         then
            current_ndist=$(( current_ndist + ndist_step ))
            echo "  $tries runs : all crashed : skipping this number of distributions"
            tries=0
            (( ncrash+=1 ))
         else
            echo "  all $local_n_chains chains crashed : retrying with more chains."
         fi
      fi
   done

	# test if success
	if [[ "$c" != "0" ]]
	then
		# store previous BIC
		previous_BIC_tilde=$BIC_tilde
      if [[ $(( $ndist % 2 )) -eq 0 ]] ; then
			previous_evenBIC_tilde=$BIC_tilde
      else
			previous_oddBIC_tilde=$BIC_tilde
      fi

		# get BIC
	   BIC_tilde=`awk NR==1 ${here}bic/dist_$ndist/bic_final_chain_$c.txt`

		echo "$ndist,$BIC_tilde" >> ${here}bic/recap_step_1.csv

      # log
      if [ $tries -eq 1 ]
      then
         echo -en " $tries run"
      else
         echo -en " $tries runs"
      fi
      echo -en " : bic=$BIC_tilde"

		# get parameters
		cp ${here}parameters/dist_$ndist/parameter_chain_$c.txt ${here}parameters/step_1/parameters_$ndist.txt

	   # we have a better BIC
		if [[ $anysuccess == false ]] || [ $ndist -eq $ndist_init ] || [[ `echo "$BIC_tilde < $BIC_star " | bc 2> /dev/null ` -eq 1 ]] ; then
		   BIC_star=$BIC_tilde
		   i=$ndist
		   nhigher=0
		   nraised=0
		   nraised_odd=0
		   nraised_even=0
		   echo -e "  : better bic found"
		else
			# increase nhigher
		   (( nhigher+=1 ))

		   # increase nraised if needed
		   if [[ `echo "$BIC_tilde > $previous_BIC_tilde " | bc 2> /dev/null ` -eq 1 ]] ; then
		      (( nraised+=1 ))
		   else
		      nraised=0
		   fi

		   # even number of dist : increase nraised if needed
		   if [[ $(( $ndist % 2 )) -eq 0 ]] ; then
		   	if [[ `echo "$BIC_tilde > $previous_evenBIC_tilde " | bc 2> /dev/null ` -eq 1 ]] ; then
			      (( nraised_even +=1 ))
			   else
			      nraised_even=0
			   fi
		   fi

		   # odd number of dist : increase nraised if needed
		   if [[ $(( $ndist % 2 )) -eq 1 ]] ; then
		   	if [[ `echo "$BIC_tilde > $previous_oddBIC_tilde " | bc 2> /dev/null ` -eq 1 ]] ; then
				   (( nraised_odd +=1 ))
				else
				   nraised_odd=0
				fi
			fi

			# log
   		echo -en "   : nhigher=$nhigher   : nraised=$nraised"
   		if [[ $transform == "mirror" ]] && [[ $mirror_bic_correction == true ]] ; then
   			echo -en "   : nraised_odd=$nraised_odd  : nraised_even=$nraised_even"
   		fi
   		echo -en "\n"
		fi


		# stopping criteria
	   if [ $ndist -ge $step1_max_dist ] ; then
	      echo -e "\nStopping criteria met : ndist >= step1_max_dist"
	      break
	   fi
	   if [ $nhigher -ge $step1_exploration_higher_limit ] ; then
	      echo -e "\nStopping criteria met : nhigher >= step1_exploration_higher_limit"
	      break
	   fi
		if [ $nraised -ge $step1_exploration_raised_limit ] ; then
		   echo -e "\nStopping criteria met : nraised >= step1_exploration_raised_limit"
		   break
		fi
	   if [[ $transform == "mirror" ]] && [[ $mirror_bic_correction == true ]] ; then
			if [ $nraised_even -ge $step1_exploration_raised_limit ] ; then
			   echo -e "\nStopping criteria met : nraised_even >= step1_exploration_raised_limit"
			   break
			fi
			if [ $nraised_odd -ge $step1_exploration_raised_limit ] ; then
			   echo -e "\nStopping criteria met : nraised_odd >= step1_exploration_raised_limit"
			   break
			fi
		fi

		# increment number of distributions
		ndist=$(( ndist + ndist_step ))

		# we have a success
		anysuccess=true

	else
		echo "all models crashed - trying no further"
		break
	fi
done


# get parameters : we will call them optimal_parameters until the end of the program
cp ${here}parameters/step_1/parameters_$i.txt ${here}parameters/optimal_parameters.txt

echo -e "\nBest BIC : $BIC_star observed with $i distributions"



###################################################################################################

if [[ "$enable_step2" == "yes" ]]
then
   #============================================
   #  STEP 2
   #============================================

   echo -e "\n----------------------------------------------------------"
   echo -e "STEP 2"
   echo -e "----------------------------------------------------------"

   # log parameters
   echo -ne "\nParameters"
   echo -ne "\n  data tranform            : $transform"
	if [[ $transform == "mirror" ]] ; then
		echo -ne "\n     Mirror parameters:"
		echo -ne "\n       symetric init    : $mirror_symetry"
		echo -ne "\n       initial nb dist  : $mirror_ndist_init"
		echo -ne "\n       step for nb dist : $mirror_ndist_step"
		echo -ne "\n       BIC correction   : $mirror_bic_correction"
	fi
   echo -ne "\n  number of resamplings    : $nb_line"
   if [[ $transform == "mirror" ]] ; then
      echo -ne " (including mirrored resamplings)"
   fi
   echo -ne "\n  model                    : $model"
   echo -ne "\n  number of chains         : $step2_chain_mode , minimal number=$step2_n_chains"
	echo -ne "\n  exploration raised limit : $step2_exploration_raised_limit"
	echo -ne "\n  exploration higher limit : $step2_exploration_higher_limit"
   echo -ne "\n  max number of crashes    : $max_EM_error"
   echo -ne "\n  EM epsilon               : $EM_epsilon"
	echo -ne "\n  EM nb iteration          : $EM_nbiteration"
	echo -ne "\n  doPreRun                 : $doPreRun"
	echo -ne "\n  PreRunInitType           : $PreRunInitType"
	echo -ne "\n  RunInitType              : $RunInitType"
   echo -ne "\n"

   # LOG
   echo -e "\n======== Iteration : 0 ========"

   echo -e "\nRunning $i distributions ..."


   # call the launch_mixture function
   # the best chain is called c
   # c is equal to zero if all chains crashed
   c=`launch_mixture $i $step2_n_chains $step2_chain_mode`

   # test if at least one chain succeeded
   if (( $c != 0 ))
   then
	   # get the bic for this best chain and i gaussians
	   # BIC_star=`awk NR==7 ${here}bic/dist_$i/bic_chain_$c.txt `
	   BIC_star=`awk NR==1 ${here}bic/dist_$i/bic_final_chain_$c.txt `

	   # comparison of BICs : step 1) and step 2)
	   if [[ `echo "$BIC_star < $BIC_tilde" | bc 2> /dev/null ` -eq 1 ]]
	   then
		   # update optimal parameters : parameters from the best of these 6 chains
		   cp ${here}parameters/dist_$i/parameter_chain_$c.txt ${here}parameters/optimal_parameters.txt
	   else
		   # get BIC : the same as in step 1)
		   BIC_star=$BIC_tilde
	   fi
   else
	   # All chains have crashed : keep the previous BIC

	   # log
	   echo -e "\nall chains have crashed"

	   # get BIC : the same as in step 1)
	   BIC_star=$BIC_tilde
   fi

   echo -e "\nOptimal BIC now : $BIC_star"


   #=================================
   #  EXPLORATION AROUND i
   #=================================

   # boolean
   go_right=true
   go_left=true

   # counters
   step_right=0
   step_left=0
   i_left=$i
   i_right=$i
   raise_left=0
   raise_right=0
   raise_even_left=0
   raise_even_right=0
   raise_odd_left=0
   raise_odd_right=0

   BIC_left=$BIC_star  # BIC on the left
   BIC_right=$BIC_star # BIC on the right
   BIC_even_left=$BIC_star  # BIC on the left
   BIC_even_right=$BIC_star # BIC on the right
   BIC_odd_left=$BIC_star  # BIC on the left
   BIC_odd_right=$BIC_star # BIC on the right
   iteration=0 	    # count the iterations

   # Save all results of this exploration on the left & right in a recap.csv file
   echo "iteration,n_distributions,gaussian_left,gaussian_right,current_BIC,\
   BIC_left,BIC_right" > ${here}bic/recap_step_3.csv

   #
   # Loop over models
   #
   while [ $go_left == true ] || [ $go_right == true ]
   do

      # increment the iteration variable
      (( iteration+=1 ))

      # LOG
      echo -e "\n======== Iteration : $iteration ========"

	   ### exploration on the left
	   if [ $go_left == true ]
	   then

         # LOG
         echo -en "\nLeft side :"

		   # increment cursor on the left
		   step_left=$(( step_left + ndist_step ))

		   # Test if the number of gaussians is larger than 1
		   g_left=$(( $i_left - $step_left ))
		   if (( $g_left >= 1 ))
		   then
	         # LOG
	         echo -e " Running ${g_left} distributions ... "

			   # Launch several chains with (i_left-step_left) gaussians
			   c_left=`launch_mixture $g_left $step2_n_chains $step2_chain_mode`

			   # check if at least 1 chain succeed
			   if (( $c_left != 0 ))
			   then
				   # store BIC_left
				   BIC_left_old=$BIC_left
					if [[ $(( $i % 2 )) -eq 0 ]] ; then
						BIC_even_left_old=$BIC_left
					else
						BIC_odd_left_old=$BIC_left
					fi

				   # get BIC_left
	            BIC_left=`awk NR==1 ${here}bic/dist_${g_left}/bic_final_chain_${c_left}.txt`

				   echo -en "   - new BIC = $BIC_left"

				   # Is it raising ?
               if [[ `echo "$BIC_left > $BIC_left_old" | bc 2> /dev/null ` -eq 1 ]] ; then
                  echo -ne "   is higher"
                  (( raise_left+=1 ))
               else
                  echo -ne "   is lower"
                  raise_left=0
               fi

					# even number of dist : Is it raising ?
					if [[ $(( $i % 2 )) -eq 0 ]] ; then
						if [[ `echo "$BIC_even_left > $BIC_even_left_old" | bc 2> /dev/null ` -eq 1 ]] ; then
							(( raise_even_left +=1 ))
						else
							raise_odd_left=0
						fi
					fi

					# odd number of dist : Is it raising ?
					if [[ $(( $i % 2 )) -eq 1 ]] ; then
						if [[ `echo "$BIC_odd_left > $BIC_odd_left_old" | bc 2> /dev/null ` -eq 1 ]] ; then
							(( raise_odd_left +=1 ))
						else
							raise_odd_left=0
						fi
					fi

					# log
               echo -ne " : raise_left now = $raise_left"
					if [[ $transform == "mirror" ]] && [[ $mirror_bic_correction == true ]] ; then
						echo -en "   : raise_odd_left=$raise_odd_left  : raise_odd_left=$raise_odd_left"
					fi
					echo -en "\n"


				   # get parameters_left
				   rm -fr ${here}parameters/parameters_step_3_left.txt
				   cp ${here}parameters/dist_$g_left/parameter_chain_$c_left.txt ${here}parameters/parameters_step_3_left.txt
			   fi
		   else
		      # g_left=0 : left limit is reached
		      echo -e " Left limit is reached"
		      go_left=false
		   fi
	   fi

	   ### exploration on the right
	   if [ $go_right == true ]
	   then

         # LOG
         echo -en "\nRight side :"

		   # increment cursor on the right
		   step_right=$(( step_right + ndist_step ))

		   # test if the number of gaussians is not too large
		   g_right=$(( $i_right + $step_right ))
		   if (( $g_right <= $max_nbr_distributions ))
		   then
	         # LOG
	         echo -e " Running ${g_right} distributions ... "

			   # launch chains with (i_right+step_right) gaussians
			   c_right=`launch_mixture $g_right $step2_n_chains $step2_chain_mode`
			   if (( $c_right != 0 ))
			   then
				   # store BIC_right
				   BIC_right_old=$BIC_right
					if [[ $(( $i % 2 )) -eq 0 ]] ; then
						BIC_even_right_old=$BIC_right
					else
						BIC_odd_right_old=$BIC_right
					fi

				   # get BIC_right
	            BIC_right=`awk NR==1 ${here}bic/dist_${g_right}/bic_final_chain_${c_right}.txt`

				   echo -ne "   - new BIC = $BIC_right"

				   # Is it raising ?
               if [[ `echo "$BIC_right > $BIC_right_old" | bc 2> /dev/null ` -eq 1 ]]
               then
                  echo -ne "   is higher"
                  (( raise_right+=1 ))
               else
                  echo -ne "   is lower"
                  raise_right=0
               fi

					# even number of dist : Is it raising ?
					if [[ $(( $i % 2 )) -eq 0 ]] ; then
						if [[ `echo "$BIC_even_right > $BIC_even_right_old" | bc 2> /dev/null ` -eq 1 ]] ; then
							(( raise_even_right +=1 ))
						else
							raise_even_right=0
						fi
					fi

					# odd number of dist : Is it raising ?
					if [[ $(( $i % 2 )) -eq 1 ]] ; then
						if [[ `echo "$BIC_odd_right > $BIC_odd_right_old" | bc 2> /dev/null ` -eq 1 ]] ; then
							(( raise_odd_right +=1 ))
						else
							raise_odd_right=0
						fi
					fi

					# log
               echo -ne " : raise_right now = $raise_right"
					if [[ $transform == "mirror" ]] && [[ $mirror_bic_correction == true ]] ; then
						echo -en "   : raise_odd_right=$raise_odd_right  : raise_odd_right=$raise_odd_right"
					fi
					echo -en "\n"

				   # get parameters_right
				   rm -fr ${here}parameters/parameters_step_3_right.txt
				   cp ${here}parameters/dist_$g_right/parameter_chain_$c_right.txt ${here}parameters/parameters_step_3_right.txt
			   fi
		   else
		      # right limit max_nbr_distributions is reached
		      echo -e " Right limit is reached"
		      go_right=false
		   fi
	   fi


	   #-----------------------
	   # Test all possibilities
	   #-----------------------

      vecho -e "\ni=$i  i_right=$i_right BIC_right=$BIC_right  i_left=$i_left BIC_left=$BIC_left"

	   # 1st possibility : BIC_left < BIC_star <= BIC_right
	   # better model on the left : stop trying on the right
	   if [[ `echo "$BIC_left < $BIC_star" | bc 2> /dev/null ` -eq 1 ]] && [[ `echo "$BIC_right >= $BIC_star" | bc 2> /dev/null` -eq 1 ]]
	   then
         # log
         vecho "1st possibility : BIC_left < BIC_star <= BIC_right"

		   # update i on the left
		   i_left=$(( i_left - step_left ))

		   # update i
		   i=$i_left

		   # update BIC
		   BIC_star=$BIC_left

		   # update the cursor on left
		   step_left=0

		   # update parameters
		   rm -fr ${here}parameters/optimal_parameters.txt
		   cp ${here}parameters/parameters_step_3_left.txt ${here}parameters/optimal_parameters.txt
	   fi

	   # 2nd possibility : BIC_right < BIC_star <= BIC_left
	   # better model on the right
	   if [[ `echo "$BIC_right < $BIC_star" | bc 2> /dev/null ` -eq 1 ]] && [[ `echo "$BIC_left >= $BIC_star" | bc 2> /dev/null ` -eq 1 ]]
	   then
         # log
         vecho "2nd possibility : BIC_right < BIC_star <= BIC_left"

		   # update i on the right
		   i_right=$(( i_right + step_right ))

		   # update i
		   i=$i_right

		   # update BIC
		   BIC_star=$BIC_right

		   # update cursor on the right
		   step_right=0

		   # update parameters for step_3
		   rm -fr ${here}parameters/optimal_parameters.txt
		   cp ${here}parameters/parameters_step_3_right.txt ${here}parameters/optimal_parameters.txt
	   fi

	   # 3rd possibility : BIC_right < BIC_star and  BIC_left < BIC_star
	   # better model on the left and right
	   if [[ `echo "$BIC_right < $BIC_star" | bc 2> /dev/null ` -eq 1 ]] && [[ `echo "$BIC_left < $BIC_star" | bc 2> /dev/null ` -eq 1 ]]
	   then
         # log
         vecho "3rd possibility : BIC_right < BIC_star and  BIC_left < BIC_star"

		   # test if BIC_left < BIC_right
		   if [ "$(echo "$BIC_left < $BIC_right" | bc)" -eq 1 ]
		   then

			   # update i on the left
			   i_left=$(( i_left - step_left ))

			   # update i
			   i=$i_left

			   # update BIC
			   BIC_star=$BIC_left

			   # update the cursor on left
			   step_left=0

			   # update parameters
			   rm -fr ${here}parameters/optimal_parameters.txt
			   cp ${here}parameters/parameters_step_3_left.txt ${here}parameters/optimal_parameters.txt
		   fi

		   # test if BIC_right < BIC_left
		   if [[ `echo "$BIC_right < $BIC_left" | bc 2> /dev/null ` -eq 1 ]]
		   then
			   # update i on the right
			   i_right=$(( i_right + step_right ))

			   # update i
			   i=$i_right

			   # update BIC
			   BIC_star=$BIC_right

			   # update cursor on the right
			   step_right=0

			   # update parameters
			   rm -fr ${here}parameters/optimal_parameters.txt
			   cp ${here}parameters/parameters_step_3_right.txt ${here}parameters/optimal_parameters.txt
		   fi

		   # test if BIC_left = BIC_right
		   if [[ `echo "$BIC_right == $BIC_left" | bc 2> /dev/null` -eq 1 ]]
		   then
			   # update i on the right and left
			   i_right=$(( i_right + step_right ))
			   i_left=$(( i_left - step_left ))

			   # update BIC on right or left
			   BIC_star=$BIC_left

			   # update cursor on the right and left
			   step_left=0
			   step_right=0

			   # update parameters
			   rm -fr ${here}parameters/optimal_parameters.txt
			   cp ${here}parameters/parameters_step_3_left.txt ${here}parameters/optimal_parameters.txt
		   fi
	   fi

      # log
      echo -e "\nCurrent best BIC : $BIC_star , observed with $i distributions"

      echo -e "   - step_left  : $step_left"
      echo -e "   - step_right : $step_right"
      echo -e "   - raise_left  : $raise_left"
      echo -e "   - raise_right : $raise_right"
      echo -e "   - go_left  : $go_left"
      echo -e "   - go_right : $go_right"
      echo -ne "\n"


      # stopping criteria
      if [ $go_left == true ] ; then
	      if (( $step_left >= $step2_exploration_higher_limit )) ; then
	         echo -e "\nStopping criteria met on the left side : step_left >= step2_exploration_higher_limit"
		      go_left=false
	      fi
	      if (( $raise_left >= $step2_exploration_raised_limit )) ; then
	         echo -e "\nStopping criteria met on the left side : raise_left >= step2_exploration_raised_limit"
		      go_left=false
	      fi
			if [[ $transform == "mirror" ]] && [[ $mirror_bic_correction == true ]] ; then
		      if (( $raise_even_left >= $step2_exploration_raised_limit )) ; then
		         echo -e "\nStopping criteria met on the left side : raise_even_left >= step2_exploration_raised_limit"
			      go_left=false
		      fi
		      if (( $raise_odd_left >= $step2_exploration_raised_limit )) ; then
		         echo -e "\nStopping criteria met on the left side : raise_odd_left >= step2_exploration_raised_limit"
			      go_left=false
		      fi
			fi
      fi
      if [ $go_right == true ] ; then
         if (( $step_right >= $step2_exploration_higher_limit )) ; then
            echo -e "\nStopping criteria met on the right side : step_right >= step2_exploration_higher_limit"
	         go_right=false
         fi
	      if (( $raise_right >= $step2_exploration_raised_limit )) ; then
	         echo -e "\nStopping criteria met on the right side : raise_right >= step2_exploration_raised_limit"
		      go_right=false
	      fi
         if [[ $transform == "mirror" ]] && [[ $mirror_bic_correction == true ]] ; then
		      if (( $raise_even_right >= $step2_exploration_raised_limit )) ; then
		         echo -e "\nStopping criteria met on the right side : raise_even_right >= step2_exploration_raised_limit"
			      go_right=false
		      fi
		      if (( $raise_odd_right >= $step2_exploration_raised_limit )) ; then
		         echo -e "\nStopping criteria met on the right side : raise_odd_right >= step2_exploration_raised_limit"
			      go_right=false
		      fi
	      fi
      fi

      # save results of exploration in a recap.csv file
      echo "$iteration,$i,$i_left,$i_right,$BIC_star,$BIC_left,$BIC_right" >> ${here}bic/recap_step_3.csv

   # end of step 2)
   done
fi

#####################################################################


#
# We now have the best gaussian model and the
# parameters associated : we can calculate the p-value
#

#----------------------------------------
# Confidence Interval for p-value
#----------------------------------------

#
# _ Use the SEM algorithm to generate confident interval
# _ Mixmod_debug_1 must be used !
# _ SEM algorithm often crashs, like the EM algorithm, since the theta
#   parameter is more or less concentrated around the local maximum of Likelihood.
#   An emplty group might appear during this algorithm
#

if [[ $doCI == "yes" ]] && [[ "$model" == "gaussian" ]] ; then

	echo -e "\n--------------------------------------"
	echo -e " SEM : pvalue confidence interval"
	echo -e "--------------------------------------\n"

	# boolean to see if we keep trying to get a confident interval
	# Variable EM_error in the options file gives us how many tries we do
	goCI=true
	tryCI=0

  # deal with null variance
   (
	cat <<- EOF
		# import parameters
		parameters = as.matrix(read.table(file="${here}parameters/optimal_parameters.txt"))

		# how many distributions ?
		K = nrow(parameters)/3
		index=1:( K * 3)

		# proportions, means & variances
		p = as.vector(parameters[which(index%%3 == 1),1])
		m = as.vector(parameters[which(index%%3 == 2),1])
		v = as.vector(parameters[which(index%%3 == 0),1])

      # deal with null variances
      v[v==0] <- 0.00000001

		# new parameters
		param2 = array(dim=length(parameters))

		j = 1
		for (i in 1:K) {
			param2[j]   = p[i]
			param2[j+1] = m[i]
			param2[j+2] = v[i]
			j = j + 3
		}

		write(param2, ncolumns=1, file="${here}parameters/optimal_parameters_nonnullv.txt")
	EOF
	) > ${jobsdir}jobs/R/null_variance.R

   ( R --no-restore --no-save --no-readline -q < ${jobsdir}jobs/R/null_variance.R ) &> ${jobsdir}jobs/R/log_null_variance_R.txt

	while [ $goCI == true ]
	do

		(
		cat <<- EOF
		NbLines
			${nb_line}
		PbDimension
			1
		NbCriterion
			1
		ListCriterion
			BIC
		NbNbCluster
			1
		ListNbCluster
			$i
		NbModel
			1
		ListModel
			Gaussian_pk_Lk_I
		NbStrategy
			1
		InitType
			USER
		InitFile
			${here}parameters/optimal_parameters_nonnullv.txt
		NbAlgorithm
			1
		Algorithm
			SEM
		StopRule
			NBITERATION
		StopRuleValue
			1100
		DataFile
			${jobsdir}inputs/mixture.dat
		EOF
		) > ${jobsdir}jobs/mixmod_sem.xem

		mixmod_debug ${jobsdir}jobs/mixmod_sem.xem > ${jobsdir}log_files/log_mixmod_sem.txt

		if ! grep "MIXMOD ERROR" ${jobsdir}log_files/log_mixmod_sem.txt >/dev/null
		then

			### Log of the SEM algorithm contains all parameters after each iteration (mode debug) ###
			mkdir ${here}parameters/SEM/ 2> /dev/null

			# get all proportions
			cat ${jobsdir}log_files/log_mixmod_sem.txt | grep "proportion" | cut -c 16-30 > ${here}parameters/SEM/proportions.txt

			# get all means
			cat ${jobsdir}log_files/log_mixmod_sem.txt | grep "mean" | cut -c 10-30 > ${here}parameters/SEM/means.txt

			# get all variances
			cat ${jobsdir}log_files/log_mixmod_sem.txt | grep "sigma" -A1 | sed -n '/sigma/!p' | sed -n '/--/!p' | cut -c 5-20 > ${here}parameters/SEM/variances.txt

			# create the R job to obtain C.I.
			(
			cat <<- EOF

			#---------------------------
			# get pvalue of maxLL SEM point
			#---------------------------

			# import parameters
			parameters=as.matrix(read.table("BICparameter.txt"))

			# how many gaussian distributions ?
			K=nrow(parameters)/3 ; index=1:(3*K) ;

			# proportions, means & variances
			proportions=as.vector(parameters[which(index%%3 == 1),1])
			proportions=proportions  + (1-sum(proportions))/length(proportions)
			means=as.vector(parameters[which(index%%3 == 2),1])
			variances=as.vector(parameters[which(index%%3 == 0),1])

         # deal with null variances
         variances[variances==0] <- 0.00000001

			# pvalue
			u = sum(proportions*pnorm(abs($quantile),means,sqrt(variances),lower.tail=FALSE)) # right
			v = sum(proportions*pnorm(-abs($quantile),means,sqrt(variances),lower.tail=TRUE)) # left

			w = sum(proportions*pnorm(1,means,sqrt(variances),lower.tail=FALSE)) # excess on the right
			z = sum(proportions*pnorm(-1,means,sqrt(variances),lower.tail=TRUE)) # excess on the left

			if ( "$transform" == "fold" ) {
				pvalue=u
				pvalue_excess_over1=w
			} else {
				pvalue=u+v
				pvalue_excess_over1=w+z
			}

			# save results
			results=matrix(0,ncol=2,nrow=1)
			results[1,1]=pvalue
			results[1,2]=pvalue_excess_over1

			colnames(results)=c("pvalue","pvalue_excess_over1")
			write.csv(results, file="${here}parameters/SEM_pvalue.csv", quote=FALSE, row.names=FALSE)

			#---------------------------
			# generation of C.I. with R
			#---------------------------

			all_proportions=as.vector(read.table('${here}parameters/SEM/proportions.txt'))
			all_means=as.vector(read.table('${here}parameters/SEM/means.txt'))
			all_variances=as.vector(read.table('${here}parameters/SEM/variances.txt'))

         # deal with null variances
         all_variances[all_variances==0] <- 0.00000001

			# get matrix
			all_proportions=as.matrix(all_proportions)
			all_means=as.matrix(all_means)
			all_variances=as.matrix(all_variances)

			# remove 100 first iterations due to the burn-in
			all_proportions=all_proportions[(102*$i+1):nrow(all_proportions),1]
			all_means=all_means[(102*$i+1):nrow(all_means),1]
			all_variances=all_variances[(102*$i+1):nrow(all_variances),1]

			# how many iterations ?
			size=length(all_proportions)/$i
			pvalues=matrix(0, ncol=1, nrow=size)

			for (j in 1:size)
			{
				# get parameters
				proportions=all_proportions[((j-1)*$i+1):((j-1)*$i+$i)]
				means=all_means[((j-1)*$i+1):((j-1)*$i+$i)]
				variances=all_variances[((j-1)*$i+1):((j-1)*$i+$i)]

				# calculate the p-value
				u = sum(proportions*pnorm($quantile,means,sqrt(variances),lower.tail=FALSE))
				v = sum(proportions*pnorm(-$quantile,means,sqrt(variances),lower.tail=TRUE))
				w = sum(proportions*pnorm(-1,means,sqrt(variances)))                 # excess on the left
				z = sum(proportions*pnorm(1,means,sqrt(variances),lower.tail=FALSE)) # excess on the right

				if ("$transform" == "fold") {
					pvalues[j]=u
				} else {
					pvalues[j]=u+v
				}

			}

			pvalues=sort(pvalues)
			q1=quantile(pvalues,0.025)
			q2=quantile(pvalues,0.975)
			pvalues=pvalues[which(pvalues>=q1)]
			pvalues=pvalues[which(pvalues<=q2)]
			IC_lower=min(pvalues)
			IC_upper=max(pvalues)
			res=matrix(0, ncol=2, nrow=1)
			res[1,1]=IC_lower
			res[1,2]=IC_upper
			colnames(res)=c('IC_lower','IC_upper')
			write.csv(res,file='${here}parameters/CI.csv', quote=FALSE, row.names=FALSE)

			EOF
			) > ${jobsdir}jobs/R/CI.R

			# R
			( R --no-restore --no-save --no-readline -q < ${jobsdir}jobs/R/CI.R ) &> ${jobsdir}jobs/R/CI_R.txt

			# success
			goCI=false
			success_CI=true

			# get SEM BIC
			SEM_bic=$(sed -n 7p BICnumericStandard.txt)
			cp BICparameter.txt ${here}parameters/SEM_parameters.txt

			# Correct BIC if mirror : half the number of independent parameters
			# BIC= -2.L + nb_free_param*ln(n)
			if [[ $transform == "mirror" ]] && [[ $mirror_bic_correction == true ]] ; then
				# ici nb_free_param=(3*K-1)
   			LogL=$(echo "scale=20; -0.5 * ( $SEMbic - (3*$i-1)*l($nb_line) ) " | bc -l )
   			if [[ $(( $i % 2 )) -eq 0 ]] ; then
   				# even
   				SEM_bic=$(echo "scale=20 ; -1 * ($LogL) + (3*$i/2-1)*l( $nb_line / 2) " | bc -l )
   			else
   				# odd
   				SEM_bic=$(echo "scale=20 ; -1 * ($LogL) + (3*($i-1)/2 + 2 - 1 )*l( $nb_line / 2) " | bc -l )
				fi
   		fi

			# add into BIC list - Chain=0 is the indicator for SEM
         echo -n "$variable,$PCA_component,$model,$i,0,$SEM_bic," >> ${here}"bics.csv"
         tail -1 ${here}parameters/SEM_pvalue.csv >> ${here}"bics.csv"

			# NO UPDATE IF BIC IS BETTER : GIVES SPURIOUS RESULTS FROM TIME TO
			# check if BIC is better
			# if [[ `echo "$SEM_bic < $BIC_star " | bc 2> /dev/null ` -eq 1 ]] ; then
			#	# SEM gives a better BIC
			#	BIC_star=$SEM_bic
			#	cp ${here}parameters/optimal_parameters.txt ${here}parameters/EM_parameters.txt
			#	cp ${here}parameters/SEM_parameters.txt ${here}parameters/optimal_parameters.txt
			#
			#	echo -e "Better BIC found : $BIC_star : parameters updated\n"
			# else
				echo -e "SEM BIC    : $SEM_bic\n"
			# fi

			# get SEM pvalue
			SEM_pvalue=`tail -1 ${here}parameters/SEM_pvalue.csv | cut --delimiter=',' --fields=1`

			# get pvalue CI
			CI_lower=`tail -1 ${here}parameters/CI.csv | cut -d',' -f1`
			CI_upper=`tail -1 ${here}parameters/CI.csv | cut -d',' -f2`

			# log
			echo -e "SEM pvalue : $SEM_pvalue"
			echo -e "        CI : $CI_lower , $CI_upper"

		else
			# increment try
			(( tryCI+=1 ))
			echo -e "$tryCI errors for generation of confidence interval"
			if [ $tryCI -gt ${max_EM_error} ]
			then
				echo -e "\nToo many errors for generation of C.I. !"

				# too many errors : stop trying
				goCI=false
				success_CI=false
				SEM_bic="/N"
				SEM_pvalue="/N"
				CI_lower="/N"
				CI_upper="/N"

			fi
		fi

	done

fi


#-------------------------------------------------------------------------
# Calculation of pvalue and graph of Mixture density
#-------------------------------------------------------------------------

echo -e "\n--------------------------------------------------------------"
echo -e "FINALIZE: Calculation of the pvalue and graph of Mixture density"
echo -e "--------------------------------------------------------------\n"

# graph file names
prefix=${here}${variable}
case "$loading_type" in
	component)
			prefix=${prefix}_comp${PCA_component}
			;;
	projection)
			prefix=${prefix}_dim${PCA_component}
			;;
esac
prefix=${prefix}_${model}
case "$transform" in
     none)
           prefix=${prefix}_raw
           ;;
     fold)
           prefix=${prefix}_folded
           ;;
     mirror)
           prefix=${prefix}_mirrored
           ;;
esac

# graph title
case "$loading_type" in
	component)
			maintitle="$title - $variable for PCA component $PCA_component"
			;;
	projection)
			maintitle="$title - $variable for $PCA_component dim projection"
			;;
esac


# R code
(
cat <<- EOF
	if ("$model" == "gaussian") {
		# import dataset
		datas=as.matrix(read.table("${jobsdir}inputs/mixture.dat"))
		abscisse=seq(-1,1,2/($nb_line-1))
		graph=array(data=c(0),dim=length(abscisse))

		# import parameters
		parameters=as.matrix(read.table("${here}parameters/optimal_parameters.txt"))

		# how many gaussian distributions ?
		K=nrow(parameters)/3 ; index=1:(3*K) ;

		# proportions, means & variances
		proportions=as.vector(parameters[which(index%%3 == 1),1])
		proportions=proportions  + (1-sum(proportions))/length(proportions)
		means=as.vector(parameters[which(index%%3 == 2),1])
		variances=as.vector(parameters[which(index%%3 == 0),1])

      # deal with null variances
      variances[variances==0] <- 0.00000001

		# pvalue
		u = sum(proportions*pnorm(abs($quantile),means,sqrt(variances),lower.tail=FALSE))
		v = sum(proportions*pnorm(-abs($quantile),means,sqrt(variances),lower.tail=TRUE))

		w = sum(proportions*pnorm(1,means,sqrt(variances),lower.tail=FALSE)) # excess on the right
		z = sum(proportions*pnorm(-1,means,sqrt(variances),lower.tail=TRUE)) # excess on the left

		if ("$transform" == "fold") {
			pvalue=u
			pvalue_excess_over1=w
		} else {
			pvalue=u+v
			pvalue_excess_over1=w+z
		}

		# graph : mixture density
		xgraph <- abscisse
		graphi=array(data=c(0),dim=c(K,length(abscisse) ) )
		for (j in 1:K)
		{
			graphi[j,] <- proportions[j]*dnorm(abscisse,means[j],sqrt(variances[j]))
			graph <- graph + graphi[j,]
		}
	}

	if ("$model" == "beta") {
		# import dataset
		datas=as.matrix(read.table("${jobsdir}inputs/mixture.dat")) ;
		abscisse=seq(0,1,1/($nb_line-1))
		graph=array(data=c(0),dim=length(abscisse))

		# import parameters
		parameters=as.matrix(read.table("${here}parameters/optimal_parameters.txt"))

		# how many beta distributions ?
		K=nrow(parameters)/3 ; index=1:(3*K) ;

		# proportions, alpha and beta
		proportions=as.vector(parameters[which(index%%3 == 1),1])
		proportions=proportions + (1-sum(proportions))/length(proportions)
		alpha=as.vector(parameters[which(index%%3 == 2),1])
		beta=as.vector(parameters[which(index%%3 == 0),1])

		if ("$transform" == "fold") {
			# pvalue
			pvalue_excess_over1=0 # by definition
			u = sum(proportions*pbeta(abs($quantile),alpha,beta,lower.tail=FALSE)) 	# right
			pvalue=u

			# graph : mixture density
			xgraph <- abscisse
			graphi=array(data=c(0),dim=c(K,length(abscisse) ) )
			for (j in 1:K)
			{
				graphi[j,] <- proportions[j]*dbeta(abscisse,alpha[j],beta[j])
				graph <- graph + graphi[j,]
			}
		} else {
			# pvalue
			pvalue_excess_over1=0 # by definition
			u = sum(proportions*pbeta(0.5+0.5*abs($quantile),alpha,beta,lower.tail=FALSE)) 	# right
			v = sum(proportions*pbeta(0.5-0.5*abs($quantile),alpha,beta,lower.tail=TRUE))		# left
			pvalue=u+v

			# graph : mixture density
			xgraph <- 2*abscisse-1
			graphi=array(data=c(0),dim=c(K,length(abscisse) ) )
			for (j in 1:K)
			{
				graphi[j,] <- proportions[j]*dbeta(abscisse,alpha[j],beta[j])/2
				graph <- graph + graphi[j,]
			}
		}
	}

	# ------------------------------------
	# save results
	results=matrix(0,ncol=5,nrow=1)
	results[1,1]=$i
	results[1,2]=$BIC_star
	results[1,3]="$model"
	results[1,4]=pvalue
	results[1,5]=pvalue_excess_over1

	colnames(results)=c("n_distibutions","BIC","model","pvalue","pvalue_excess_over1")
	write.csv(results, file="${here}results.csv", quote=FALSE, row.names=FALSE)

	if ( "$dograph" == "yes" ) {

		# ------------------------------------
		# Graph parameters
		source("$R_dir/standard_plot_functions.R")
		source("$R_dir/settings_highres.R")
		source("$R_dir/standard_settings_highdensity.R")

		cex.main=1.3

		# ------------------------------------
		# Main graph : data histogram and mixture model

		grdev.res    = 100
		grdev.width  = 1024
		grdev.height = 768
		grdev.type   = "png" # tiff, png, svg ...

		initgraphs("${prefix}_mixture.png")

		# if (abs($quantile)<1) {
		#    xlim=c(-1,1)
		# } else {
			if ("$transform" == "fold") {
				xlim=abs($quantile)*1.5
				xlim=c(0,xlim)
			} else {
				xlim=abs($quantile)*1.3
				xlim=c(-xlim,xlim)
			}
		# }

		maxgraph <- max( graph[graph != Inf] )
		ylim=c(0,maxgraph )

		if ("$loading_type" == "component") {
			xlab="Loading $PCA_component"
		} else {
			xlab="Projected Loading $PCA_component dimension"
		}

		hist(datas, freq=FALSE, xlab=xlab, ylab="Density", xlim=xlim , col="green", border="darkgreen"
				# , breaks = length(abscisse) / 100
				, breaks = "Scott"
			  , ylim=ylim, main="")

		lines(xgraph, graph , lty=1, col="darkblue", lwd=1.2)

		lines(x=c(abs($quantile),abs($quantile)) , y=c(0,maxgraph) , col="red", lwd=1, lty=2)

		if ("$transform" == "fold") {
			lines(x=c(-1*abs($quantile),-1*abs($quantile)) , y=c(0,maxgraph) , col="red", lwd=1, lty=2)
		}

		if ("$loading_type" == "component") {
			leg.txt=c("Histogram of resampled loadings","${model} mixture density",paste("Observed loading$PCA_component =",formatC($quantile,format="f",digits=2)))
		} else {
			leg.txt=c("Histogram of resampled loadings","${model} mixture density",paste("Observed proj. loading dim$PCA_component =",formatC($quantile,format="f",digits=2)))
		}
		legend("topright", legend=leg.txt, cex=0.8, col=c("green","darkblue","red"), pch="", lty=c(1,1,1), lwd=2
				  , box.lwd = 0.5, box.lty = 1, box.col = "grey50")

		#add title - PCA_component
		main="$maintitle"
		subtitle=paste("Number of $model distributions=" , K, "  BIC=", signif($BIC_star,digits=5) , "   p-value=",signif(pvalue,digits=3), sep="")
		mtext(text=main,side=3,line=2,cex=cex.main)
		mtext(text=subtitle,side=3, line=1,cex = cex.sub)

		# expand XY axes
		box(bty="l",col="black")

		closegraphs()

		# ------------------------------------
		# distributions with sample density

		initgraphs("${prefix}_distributions_density.png")

		hdata <- hist(datas, breaks = "Scott", plot=FALSE )

		plot(   hdata\$mids, hdata\$density
				 , main=""
				 , xlab = "", ylab = ""
				 , xlim=xlim
				 , ylim=ylim
				 , col="darkgreen"
				 , bg=rgb(0,0,255, 10, maxColorValue=255)
				 , pch=3 , cex=1
				 , axes = FALSE
			 )

		lines(hdata\$mids, hdata\$density, lty=1, col="darkgreen", lwd=2)

		lines(xgraph, graph , lty=1, col="darkblue", lwd=1.5)

		for (j in 1:K) {
			lines(xgraph, graphi[j,] , lty=2, col="black", lwd=0.5)
		}

		lines(x=c(abs($quantile),abs($quantile)) , y=c(0,maxgraph) , col="red", lwd=1, lty=2)

		if ("$transform" == "fold") {
			lines(x=c(-1*abs($quantile),-1*abs($quantile)) , y=c(0,maxgraph) , col="red", lwd=1, lty=2)
		}

		if ("$loading_type" == "component") {
			leg.txt=c("Resampled loadings density","${model} mixture density","Distributions",paste("Observed loading$PCA_component =",formatC($quantile,format="f",digits=2)))
		} else {
			leg.txt=c("Resampled loadings density","${model} mixture density","Distributions",paste("Observed proj. loading dim$PCA_component =",formatC($quantile,format="f",digits=2)))
		}

		legend("topright", legend=leg.txt, cex=0.8, col=c("green","darkblue","black","red"), pch="", lty=c(1,1,2,2), lwd=c(2,2,1,1)
				  , box.lwd = 0.5, box.lty = 1, box.col = "grey50")

		dist.txt <- c("model     : $model", "transform : $transform" , "size      : $nb_line", " ")

		if ("$model" == "gaussian") {
			dist.txt <- c(dist.txt,"distributions : " , "proportion, mean, variance")
			for (j in 1:K) {
				dist.txt <- c(dist.txt , paste( j , " : p="
					                             , formatC(proportions[j],format="f",digits=3)
					                             , " m="
					                             , formatC(means[j],format="f",digits=2,flag="+")
					                             , " v="
					                             , formatC(variances[j],format="g",digits=2)
					                             , sep=""
					                           )
					          )
			}
		}
		if ("$model" == "beta") {
			dist.txt <- c(dist.txt, "distributions : " , "proportion, alpha, beta")
			for (j in 1:K) {
				dist.txt <- c(dist.txt , paste( j , " : p="
					                             , formatC(proportions[j],format="f",digits=3)
					                             , " a="
					                             , formatC(alpha[j],format="f",digits=1, width=5)
					                             , " b="
					                             , formatC(beta[j],format="f",digits=1, width=5)
					                             , sep=""
					                           )
					          )
			}
		}
		savefamily <- par(family="mono")
		savefont <- par(font=2)

		if ("$transform" == "fold") {
			legend("right", legend=dist.txt, cex=0.7, col=c("black"), pch="", lty=1, lwd=0 , box.lwd = 0)
		} else {
			legend("topleft", legend=dist.txt, cex=0.7, col=c("black"), pch="", lty=1, lwd=0 , box.lwd = 0)
		}

		par(family=savefamily)
		par(font=savefont)

		#add title - PCA_component
		main="$maintitle"
		subtitle=paste("Number of $model distributions=" , K, "  BIC=", signif($BIC_star,digits=5) , "   p-value=",signif(pvalue,digits=3), sep="")
		mtext(text=main,side=3,line=2,cex=cex.main)
		mtext(text=subtitle,side=3, line=1,cex = cex.sub)

		# axes
		axis(side=1)
		if ("$loading_type" == "component") {
			mtext(text="Loading $PCA_component" , cex=cex.axis.label , side=1 , line=3 , las=0)
		} else {
			mtext(text="Proj. Loading $PCA_component dim" , cex=cex.axis.label , side=1 , line=3 , las=0)
		}
		axis(side=2)
		mtext(text="Density", cex=cex.axis.label , side=2 , line=2, las=0)

		# expand XY axes
		box(bty="l",col="black")

		closegraphs()

		# ------------------------------------
		# distributions with sample histogram
		initgraphs("${prefix}_distributions.png")

		hist(datas, freq=FALSE
				, xlab="", ylab=""
				, xlim=xlim
				, ylim=ylim
				, col=rgb(50,255,50,     100, maxColorValue=255)
				, border=rgb(100,100,100,     150, maxColorValue=255)
				# , breaks = length(abscisse) / 100
				, breaks = "Scott"
				# , breaks = "Sturges"
				, main=""
				, axes = FALSE
			  )

		lines(xgraph, graph , lty=1, col="darkblue", lwd=1.5)

		for (j in 1:K) {
			lines(xgraph, graphi[j,] , lty=2, col="black", lwd=1.0)
		}

		lines(x=c(abs($quantile),abs($quantile)) , y=c(0,maxgraph) , col="red", lwd=1, lty=2)

		if ("$transform" == "fold") {
			lines(x=c(-1*abs($quantile),-1*abs($quantile)) , y=c(0,maxgraph) , col="red", lwd=1, lty=2)
		}

		if ("$loading_type" == "component") {
			leg.txt=c("Resampled loadings density","${model} mixture density","Distributions",paste("Observed loading$PCA_component =",formatC($quantile,format="f",digits=2)))
		} else {
			leg.txt=c("Resampled loadings density","${model} mixture density","Distributions",paste("Observed proj. loading dim$PCA_component =",formatC($quantile,format="f",digits=2)))
		}

		legend("topright", legend=leg.txt, cex=0.8, col=c("green","darkblue","black","red"), pch="", lty=c(1,1,2,2), lwd=c(2,2,1,1)
				  , box.lwd = 0.5, box.lty = 1, box.col = "grey50")

		dist.txt <- c("model     : $model", "transform : $transform" , "size      : $nb_line", " ")

		if ("$model" == "gaussian") {
			dist.txt <- c(dist.txt,"distributions : " , "proportion, mean, variance")
			for (j in 1:K) {
				dist.txt <- c(dist.txt , paste( j , " : p="
					                             , formatC(proportions[j],format="f",digits=3)
					                             , " m="
					                             , formatC(means[j],format="f",digits=2,flag="+")
					                             , " v="
					                             , formatC(variances[j],format="g",digits=2)
					                             , sep=""
					                           )
					          )
			}
		}
		if ("$model" == "beta") {
			dist.txt <- c(dist.txt, "distributions : " , "proportion, alpha, beta")
			for (j in 1:K) {
				dist.txt <- c(dist.txt , paste( j , " : p="
					                             , formatC(proportions[j],format="f",digits=3)
					                             , " a="
					                             , formatC(alpha[j],format="f",digits=1, width=5)
					                             , " b="
					                             , formatC(beta[j],format="f",digits=1, width=5)
					                             , sep=""
					                           )
					          )
			}
		}
		savefamily <- par(family="mono")
		savefont <- par(font=2)

		if ("$transform" == "fold") {
			legend("right", legend=dist.txt, cex=0.7, col=c("black"), pch="", lty=1, lwd=0 , box.lwd = 0)
		} else {
			legend("topleft", legend=dist.txt, cex=0.7, col=c("black"), pch="", lty=1, lwd=0 , box.lwd = 0)
		}

		par(family=savefamily)
		par(font=savefont)

		#add title - PCA_component
		main="$maintitle"
		subtitle=paste("Number of $model distributions=" , K, "  BIC=", signif($BIC_star,digits=5) , "   p-value=",signif(pvalue,digits=3), sep="")
		mtext(text=main,side=3,line=2,cex=cex.main)
		mtext(text=subtitle,side=3, line=1,cex = cex.sub)

		# axes
		axis(side=1)
		if ("$loading_type" == "component") {
			mtext(text="Loading $PCA_component" , cex=cex.axis.label , side=1 , line=3 , las=0)
		} else {
			mtext(text="Proj. Loading $PCA_component dim" , cex=cex.axis.label , side=1 , line=3 , las=0)
		}

		axis(side=2)
		mtext(text="Density", cex=cex.axis.label , side=2 , line=2, las=0)

		# expand XY axes
		box(bty="l",col="black")

		closegraphs()

		# ------------------------------------
		# graph bics by number of distributions
		#
		bics <- read.csv(file="${here}bics.csv",header=T)

		initgraphs("${prefix}_bics.png")

		# make room depending on size of string
		offsetY = round(log10(max(abs(bics\$bic))))
		par(mar=par()\$mar+c(0,max(offsetY - par()\$mar[2],1),0,0))

		delta=0.10*abs( max(bics\$bic)-min(bics\$bic) )
		ylim=c( min(bics\$bic)-delta , max(bics\$bic)+delta )

		plot(   bics\$n_distributions, bics\$bic
				 , main=""
				 , xlab = "", ylab = ""
				 , ylim=ylim
				 , col=rgb(0,0,255, 255, maxColorValue=255)
				 , bg=rgb(0,0,255, 10, maxColorValue=255)
				 , pch=21 , cex=1.3
				 , axes = FALSE
			 )

		par(new = TRUE)
		plot(   results[1,1], results[1,2]
				 , main=""
				 , xlab = "", ylab = ""
				 , ylim=ylim , xlim=c(min(bics\$n_distributions),max(bics\$n_distributions))
				 , col="red"
		#      , bg=bg.variable
				 , pch=3 , cex=3.0
				 , axes = FALSE
			 )

		# axes
		axis(side=1, at=as.numeric(names(table(bics\$n_distributions))) )
		mtext(text="Number of $model distributions" , cex=cex.axis.label , side=1 , line=3 , las=0)
		axis(side=2)
		mtext(text="BIC", cex=cex.axis.label , side=2 , line=offsetY , las=0)


		#add title - PCA_component
		main="$maintitle"
		subtitle=paste("Best BIC=", signif($BIC_star,digits=5) , "   Number of $model distributions=" , K , sep="")
		mtext(text=main,side=3,line=2,cex=cex.main)
		mtext(text=subtitle,side=3, line=1,cex = cex.sub)


		# expand XY axes
		box(bty="l",col="black")

		closegraphs()

		# ------------------------------------------
		# graph pvalues by number of distributions
		#

		initgraphs("${prefix}_pvalues.png")

		# log pvalues
		lp <- -log10(bics\$pvalue)

		# make room on the left
		par(mar=par()\$mar+c(0,2,0,0))

		ylim=c( 0 , max(lp) + 1 )

		plot(   bics\$n_distributions, lp
				 , main=""
				 , xlab = "", ylab = ""
				 , ylim=ylim
				 , col=rgb(0,0,255, 255, maxColorValue=255)
				 , bg=rgb(0,0,255, 10, maxColorValue=255)
				 , pch=22 , cex=1.3
				 , axes = FALSE
			 )

		# best model
		par(new = TRUE)
		plot(   results[1,1], -log10( pvalue )
				 , main=""
				 , xlab = "", ylab = ""
				 , ylim=ylim , xlim=c(min(bics\$n_distributions),max(bics\$n_distributions))
				 , col="red"
		#      , bg=bg.variable
				 , pch=3 , cex=4.0
				 , axes = FALSE
			 )

		# lowest pvalue
		par(new = TRUE)
		plot(  bics[ bics\$pvalue == min(bics\$pvalue) , ]\$n_distributions[1] , -log10( min(bics\$pvalue) )
				 , main=""
				 , xlab = "", ylab = ""
				 , ylim=ylim , xlim=c(min(bics\$n_distributions),max(bics\$n_distributions))
				 , col="blue"
		#      , bg=bg.variable
				 , pch=3 , cex=3.0
				 , axes = FALSE
			 )


		# bottom axis
		axis(side=1, at=as.numeric(names(table(bics\$n_distributions))) )
		mtext(text="Number of $model distributions" , cex=cex.axis.label , side=1 , line=3 , las=0)

		# left axis
		min.p<- max(ylim,lp)
		left.axis<-c(-log10(1),-log10(0.1),-log10(0.01))
		while ( left.axis[length(left.axis)] < min.p ) {
			left.axis<-c( left.axis, left.axis[length(left.axis)]+1 )
		}

		if (length(left.axis)>10) {
			left.axis<-c(ceiling( seq(left.axis[1],left.axis[length(left.axis)], length.out=10) ) )
		}

		left.axis.label<-signif(10^-(left.axis),2)

		axis(tck=-0.025,side=2, at=left.axis, labels=left.axis.label, las=2)

		mtext(text="pvalue (log scale)", cex=cex.axis.label , side=2 , line=5 , las=0)


		#add title - PCA_component
		main="$maintitle"
		subtitle=paste("Red cross represents best model : number of $model distributions=" , K, " , p-value=",signif(pvalue,digits=3),   sep="")
		subtitle2=paste("( Blue cross represents lowest p-value=",signif(min(bics\$pvalue),digits=3), " )",  sep="")
		mtext(text=main,side=3,line=2,cex=cex.main)
		mtext(text=subtitle,side=3, line=1,cex = cex.sub)
		mtext(text=subtitle2,side=3, line=0,cex = cex.sub)


		# expand XY axes
		box(bty="l",col="black")

		closegraphs()


		# ------------------------------------------
		# graph pvalues by bics
		#

		initgraphs("${prefix}_pvalues_by_bics.png")

		# log pvalues
		lp <- -log10(bics\$pvalue)

		# make room on the left
		par(mar=par()\$mar+c(0,2,0,0))

		ylim=c( 0 , max(lp) + 1 )

		delta=0.10*abs( max(bics\$bic)-min(bics\$bic) )
		xlim=c( min(bics\$bic)-delta , max(bics\$bic)+delta )

		plot(   bics\$bic, lp
				 , main=""
				 , xlab = "", ylab = ""
				 , ylim=ylim, xlim=xlim
				 , col=rgb(0,0,255, 255, maxColorValue=255)
				 , bg=rgb(0,0,255, 10, maxColorValue=255)
				 , pch=22 , cex=1.3
				 , axes = FALSE
			 )

		# best model
		par(new = TRUE)
		plot(   results[1,2], -log10( pvalue )
				 , main=""
				 , xlab = "", ylab = ""
				 , ylim=ylim , xlim=xlim
				 , col="red"
		#      , bg=bg.variable
				 , pch=3 , cex=4.0
				 , axes = FALSE
			 )

		# lowest pvalue
		par(new = TRUE)
		plot(  bics[ bics\$pvalue == min(bics\$pvalue) , ]\$bic[1] , -log10( min(bics\$pvalue) )
				 , main=""
				 , xlab = "", ylab = ""
				 , ylim=ylim , xlim=xlim
				 , col="blue"
		#      , bg=bg.variable
				 , pch=3 , cex=3.0
				 , axes = FALSE
			 )


		# bottom axis
		axis(side=1 )
		mtext(text="BICs" , cex=cex.axis.label , side=1 , line=3 , las=0)

		# left axis
		min.p<- max(ylim,lp)
		left.axis<-c(-log10(1),-log10(0.1),-log10(0.01))
		while ( left.axis[length(left.axis)] < min.p ) {
			left.axis<-c( left.axis, left.axis[length(left.axis)]+1 )
		}

		if (length(left.axis)>10) {
			left.axis<-c(ceiling( seq(left.axis[1],left.axis[length(left.axis)], length.out=10) ) )
		}

		left.axis.label<-signif(10^-(left.axis),2)

		axis(tck=-0.025,side=2, at=left.axis, labels=left.axis.label, las=2)

		mtext(text="pvalue (log scale)", cex=cex.axis.label , side=2 , line=5 , las=0)


		#add title - PCA_component
		main="$maintitle"
		subtitle=paste("Red cross represents best model : number of $model distributions=" , K, " , BIC=", signif($BIC_star,digits=5) , " , p-value=",signif(pvalue,digits=3),  sep="")
		subtitle2=paste("( Blue cross represents lowest p-value=",signif(min(bics\$pvalue),digits=3), " )",  sep="")
		mtext(text=main,side=3,line=2,cex=cex.main)
		mtext(text=subtitle,side=3, line=1,cex = cex.sub)
		mtext(text=subtitle2,side=3, line=0,cex = cex.sub)


		# expand XY axes
		box(bty="l",col="black")

		closegraphs()

	} # if dograph


EOF
) > ${jobsdir}jobs/R/pvalue_and_graph.R

echo -e "\nP-value and graph computation are running...\n"

# R
( R --no-restore --no-save --no-readline -q < ${jobsdir}jobs/R/pvalue_and_graph.R ) &> ${jobsdir}jobs/R/log_pvalue_and_graph_R.txt

# REMOVE ALL THESE FILES (WE DON'T NEED THEM ANYMORE)
rm -fr CV* DCV* ICL* NEC* numericComplete.txt complete.txt error* BIC*
rm -fr ${here}parameters/parameters_step_3_left.txt ${here}parameters/parameters_step_3_right.txt

# TOTAL ELAPSED TIME
end=$(date +%s)
tet=$(( $end - $begin ))

echo -e "\n========================================================"
echo -e "TOTAL ELAPSED TIME :  $tet   seconds"
echo -e "========================================================\n"

##### RESULTS #####

# write result in a file
pvalue=`tail -1 ${here}results.csv | cut -d',' -f4`
pvalue_excess=`tail -1 ${here}results.csv | cut -d',' -f5`
echo "variable,component,n_distributions,BIC,pvalue,pvalue_excess_over1,SEM_BIC,SEM_pvalue,SEM_pvalue_lowerCI,SEM_pvalue_upperCI,duration" > ${here}results.csv
echo "$variable,$PCA_component,$i,$BIC_star,$pvalue,$pvalue_excess,$SEM_bic,$SEM_pvalue,$CI_lower,$CI_upper,$tet" >> ${here}results.csv

# echo final results
echo -e "\n================================================\n"
echo -e "	 +++ FINAL RESULTS +++ \n"
echo -e " Title         : $title"
echo -e " Variable      : $variable"
echo -e " PCA_component : $PCA_component"
echo -e "\n Model         : $model"
echo -e "\n Number of distributions   : $i"
echo -e " BIC criterion             : $BIC_star"
echo -e " pvalue                    : $pvalue"
echo -e " pvalue 95%CI              : [ $CI_lower ; $CI_upper ]"
echo -e "\n================================================\n"


) 2>&1 | tee ${finaljobsdir}log_all.txt

# store results (from ramdisk)
if [[ "$minimize_IO" == "yes" ]]
then
	# store only what we need
	cp $here/results.csv $finalhere
	cp $here/bics.csv $finalhere
	cp $here/*.png $finalhere
	mkdir -p $finalhere/parameters/
	cp $here/parameters/optimal_parameters.txt $finalhere/parameters/optimal_parameters.txt
else
	# store all working area
	cp -R $jobsdir/* $finaljobsdir
	# store all results
	cp -R $here/* $finalhere
fi

# clean ramdisk
rm -rf /ramdisk${finaljobsdir} # $jobsdir  # workdir
rm -rf /ramdisk${finalhere}    # $here	  # result dir

### END OF PROGRAM ###
