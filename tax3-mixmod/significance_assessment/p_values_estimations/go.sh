#!/bin/bash

####################################################
#
# launch batches for selected variables and components
#                   
# 1 parameter:     
#  $1 = options file
#                   
####################################################

echo -e "\n ---- go script ----\n"

# --------------------------------------------------------------
# check parameters
#
usage () {
  echo "Usage:"
  echo "go.sh - launch EM algorithm with gaussian or beta distribution"
  echo "syntax:"
  echo "    go.sh options_file"
  exit
} 
if [ $# -ne 1 ] 
then    
	usage
fi

# --------------------------------------------------------------
# get analysis options

# check whether options file is readable
if [ ! -r "$1" ] ; then echo "Error: options file is not readable"; exit; fi

# source the option file
. $1

# main analysis mode
model=$MODEL
case "$model" in
   beta)
      runner=runner_mixture.sh; beta="yes" ;;
   gaussian)
      runner=runner_mixture.sh; beta="no" ;;
   gaussianCI)
      runner=runner_gaussian_mixture_debug_1.sh; beta="no" ;;
   *)
      echo "option file : wrong MODEL"
      exit
esac

# check whether runner has been set
if [ -z "$runner" ] ; then echo "Error: missing EM algorithm distribution"; usage; fi;

# EM options
EM_epsilon=$EM_EPSILON
EM_nbiteration=$EM_NBITERATION

# launch options
exec_mode=$EXEC_MODE
qsub_args=$QSUB_ARGS
qsub_batch=$QSUB_BATCH
user=$(id -nu)

# get directory of analysis
here=$HERE
if expr "$here" : ".*\(.\)" != '/'>/dev/null ; then here=$here"/"; fi

# get loadings.csv
loadings=$LOADINGS

# get their type : component or projection
loading_type=$LOADING_TYPE

# get resampling results directory
results=$RESULTS
if expr "$results" : ".*\(.\)" != '/'>/dev/null ; then results=$results"/"; fi

# get output directory
outdir=$OUTDIR
if expr "$outdir" : ".*\(.\)" != '/'>/dev/null ; then outdir=$outdir"/"; fi

# get Batch directory
batchdir=$BATCH_DIR
if expr "$batchdir" : ".*\(.\)" != '/'>/dev/null ; then batchdir=$batchdir"/"; fi

# get R directory
Rdir=$R_DIR
if expr "$Rdir" : ".*\(.\)" != '/'>/dev/null ; then Rdir=$Rdir"/"; fi

# get list of variables
variables=$VARIABLES
# check whether variables file is readable
if [ ! -r "$variables" ] ; then echo "Error: variables file is not readable"; exit; fi
# create file name to store list of unique variables
var_list=$(basename variables | cut -f 1 -d '.')"_uniq.txt"

# get number of resamplings to keep for mixmod modeling
n_resamplings=$N_RESAMPLINGS

# get number of EM error
em_err=$EM_ERROR

# fixed or adaptive number of chains
step1_chain_mode=$STEP1_CHAIN_MODE
step2_chain_mode=$STEP2_CHAIN_MODE

# get number of independent chains to run
step1_n_chains=$STEP1_N_CHAINS
step2_n_chains=$STEP2_N_CHAINS

# step1 exploration
step1_exploration_raised_limit=$STEP1_EXPLORATION_RAISED_LIMIT
step1_exploration_higher_limit=$STEP1_EXPLORATION_HIGHER_LIMIT
step1_max_dist=$STEP1_MAX_DIST

# step 2
enable_step2=$ENABLE_STEP2

# step2 exploration
step2_exploration_raised_limit=$STEP2_EXPLORATION_RAISED_LIMIT
step2_exploration_higher_limit=$STEP2_EXPLORATION_HIGHER_LIMIT

# are we doing graphs? (no: speed-up?)
dograph=$DOGRAPH

# title for graphs
title=$TITLE

# multi thread
RunChainsInParallel=$RUN_CHAINS_IN_PARALLEL

#
# Transform options
#

# none	: no data transform
# mirror : duplicate data with oposite values
# fold	: take absolute values
transform=$TRANSFORM


mirror=$MIRROR					# mirror dataset
mirror_symetry=$MIRROR_SYMETRY		# make init parameters symetric
mirror_ndist_init=$MIRROR_NDIST_INIT		# number of distributions to start from
mirror_ndist_step=$MIRROR_NDIST_STEP		# step when increasing number of distributions
mirror_bic_correction=$MIRROR_BIC_CORRECTION	# BIC correction



# ----------------------------------------


# ----------------------------------------
# cut the header from list of variables
# make sure there is no (var,comp) duplicate
cat $variables | sed '1d' | sort | uniq > ${here}$var_list

# ----------------------------------------
# main loop
echo -e "\n Launching analyses ... \n"

# create directory where all results and working files are located
mkdir -p $outdir

# file listing all analyses and location (used in data management)
case "${loading_type}" in
	component)
		echo "variable,component,jobid,path" > ${outdir}analyses.csv
	   ;;
	projection)
		echo "variable,dimension,jobid,path" > ${outdir}analyses.csv	   
		;;
esac

# file for qsub job removal
rm ${here}qsub_del_all.sh 2> /dev/null
touch ${here}qsub_del_all.sh
chmod u+x ${here}qsub_del_all.sh

# temp log file
rm ${here}_temp.txt 2> /dev/null
touch ${here}_temp.txt

if [ "$qsub_batch" != "" ]
then
	# qsub command go to a batch file
	(
	cat <<- EOF
		#!/bin/bash
		
		echo -e "\n Launching analyses ... \n"

		# create directory where all results and working files are located
		mkdir -p $outdir

		# file listing all analyses and location (used in data management)
		case "${loading_type}" in
			component)
				echo "variable,component,jobid,path" > ${outdir}analyses.csv
				;;
			projection)
				echo "variable,dimension,jobid,path" > ${outdir}analyses.csv	   
				;;
		esac

		# file for qsub job removal
		rm ${here}qsub_del_all.sh 2> /dev/null
		touch ${here}qsub_del_all.sh
		chmod u+x ${here}qsub_del_all.sh

		# temp log file
		rm ${here}_temp.txt 2> /dev/null
		touch ${here}_temp.txt
		
	EOF
	) > $qsub_batch
	chmod u+x $qsub_batch
fi

# for each (var,comp)
# FORK THIS ???
ijob=0
cat ${here}$var_list | while read varcomp ; 
do
	ijob=$((1+$ijob))
	
   # get varname and component to analyse
   var=`echo $varcomp | cut --delimiter=, --fields=1 | sed -e"s/ //g"`
   comp=`echo $varcomp | cut --delimiter=, --fields=2 | sed -e"s/ //g"`
 
   # where to store results
   case "${loading_type}" in
		component)
		   vardir="comp"$comp"_"$var"/" 
			# log
			echo -ne "${ijob}  "
			printf "%15s" $var
			budimcomp='", component " $comp "   ...   "'
			echo -ne ", component " $comp "   ...   "		   
		   ;;
		projection)
		   vardir="dim"$comp"_"$var"/" 
			# log
			echo -ne "${ijob}  "
			printf "%15s" $var
			budimcomp='", dimension " $comp "   ...   "'
			echo -ne ", dimension " $comp "   ...   "		   
		   ;;
	esac

	# make working directory
	resdir=$outdir$vardir"results"
	workdir=$outdir$vardir"workdir"
	# first remove existing directory
	rm -fr $workdir $resdir
	# then create working directory
	mkdir -p $workdir $resdir

	# copy header and observed loading
   # csv file, spaces allowed
   # "rs1234 ,0.5"
	head -1 $loadings    | cut --delimiter=, --fields=1,$((1+$comp)) >  $resdir/loading.csv
	grep "^${var}\( \)\?\," $loadings  | cut --delimiter=, --fields=1,$((1+$comp)) >> $resdir/loading.csv

   # check we have only 2 lines:
   if [[ `wc -l $resdir/loading.csv | cut --fields=1 --delimiter=' '` != '2' ]] ; then
      echo "error: loading file does not have 2 lines"
      cat $resdir/loading.csv
      exit 1
   fi

	# copy resamplings
	n=$((1+$n_resamplings))
   # compcol=$(( 2 + ($comp-1)*4 ))  # we have to deal with the p & CI columns
   compcol=$(( $comp + 1 ))  # no p & CI in C++
	head -$n ${results}${var}.csv | cut --delimiter=, --fields=1,$compcol >  $resdir/$var.csv

	#
	# copy and personalize runner batch
	#
	# have first to escape slashes for sed
	escpathbatch=$(printf "$batchdir" | sed 's/\//\\\//g')
	escpathwork=$(printf "$workdir" | sed 's/\//\\\//g')
	escpathres=$(printf "$resdir" | sed 's/\//\\\//g')
	escpathR=$(printf "$Rdir" | sed 's/\//\\\//g')
	# add trailing slash
	escpathwork=$escpathwork"\/"
	escpathres=$escpathres"\/"

	# copy runner batch and update its options
   sed 	-e"s/<MODEL>/"$model"/"                   \
			\
			-e"s/<EM_EPSILON>/"$EM_epsilon"/"                   \
			-e"s/<EM_NBITERATION>/"$EM_nbiteration"/"                   \
			\
			-e"s/<VAR>/"$var"/"                   \
			-e"s/<LOADING_TYPE>/"$loading_type"/"                   \
			\
			-e"s/<HERE>/"$escpathres"/"           \
			-e"s/<JOBSDIR>/"$escpathwork"/"       \
			-e"s/<BATCH_DIR>/"$escpathbatch"/"       \
			-e"s/<RDIR>/"$escpathR"/"             \
			\
			-e"s/<STEP1_CHAIN_MODE>/"$step1_chain_mode"/"     \
			-e"s/<STEP1_N_CHAINS>/"$step1_n_chains"/"            \
			-e"s/<STEP1_EXPLORATION_RAISED_LIMIT>/"$step1_exploration_raised_limit"/"   \
			-e"s/<STEP1_EXPLORATION_HIGHER_LIMIT>/"$step1_exploration_higher_limit"/"   \
			-e"s/<STEP1_MAX_DIST>/"$step1_max_dist"/"   \
			\
			-e"s/<ENABLE_STEP2>/"$enable_step2"/"   \
			-e"s/<STEP2_CHAIN_MODE>/"$step2_chain_mode"/"     \
			-e"s/<STEP2_N_CHAINS>/"$step2_n_chains"/"            \
			-e"s/<STEP2_EXPLORATION_RAISED_LIMIT>/"$step2_exploration_raised_limit"/"   \
			-e"s/<STEP2_EXPLORATION_HIGHER_LIMIT>/"$step2_exploration_higher_limit"/"   \
			\
			-e"s/<TRANSFORM>/"$transform"/"   \
			\
			-e"s/<MIRROR_SYMETRY>/"$mirror_symetry"/"   \
			-e"s/<MIRROR_NDIST_INIT>/"$mirror_ndist_init"/"   \
			-e"s/<MIRROR_NDIST_STEP>/"$mirror_ndist_step"/"   \
			-e"s/<MIRROR_BIC_CORRECTION>/"$mirror_bic_correction"/"   \
			\
			-e"s/<EM_ERROR>/"$em_err"/"           \
			-e"s/<DOGRAPH>/$dograph/g"                 \
			-e"s/<RUN_CHAINS_IN_PARALLEL>/${RunChainsInParallel}/"                 \
			-e"s/<TITLE>/$title/g"                 \
			-e"s/<PCA_COMPONENT>/"$comp"/"        \
			\
			${batchdir}$runner > $workdir/$runner

			
	chmod u+x $workdir/$runner

	# EM algorithm needed if beta model
	cp ${batchdir}EM_functions.R $workdir; chmod u+x $workdir/EM_functions.R

	# launch batch job 
	if [ "$exec_mode" = "qsub" ]
	then 
		if [ "$qsub_batch" = "" ]
		then
			# qsub
			qsub $qsub_args -o $workdir -e $workdir $workdir/$runner | tee ${here}_temp.txt
		   echo -n "qdel " >> ${here}qsub_del_all.sh
		   cat ${here}_temp.txt >> ${here}qsub_del_all.sh
		   echo >> ${here}qsub_del_all.sh
			# log   
			jobid=`head -1 ${here}_temp.txt`
			echo $var","$comp","$jobid","$resdir >> ${outdir}analyses.csv
	   else
			# qsub command go to a batch file
			(
			cat <<- EOF
				# submit new job
				echo -ne "${ijob}  "
				printf "%15s" $var		
				echo -ne $budimcomp		
				qsub $qsub_args -o $workdir -e $workdir $workdir/$runner | tee ${here}_temp.txt
				echo -n "qdel " >> ${here}qsub_del_all.sh
				cat ${here}_temp.txt >> ${here}qsub_del_all.sh
				echo >> ${here}qsub_del_all.sh	   
				# log   
				jobid=\`head -1 ${here}_temp.txt\`
				echo $var","$comp","\$jobid","$resdir >> ${outdir}analyses.csv
				echo -ne $jobid
							
			EOF
			) >> $qsub_batch
			echo " -> batch"
	   fi
	else
		# at command
		at -f $workdir/$runner now  | tee ${here}_temp.txt
		# log   
		jobid=`head -1 ${here}_temp.txt`
		echo $var","$comp","$jobid","$resdir >> ${outdir}analyses.csv
	fi


done
 

if [ "$exec_mode" = "qsub" ] && [ "$qsub_batch" != "" ]
then
	(
	cat <<- EOF
		# clean up
		rm -f ${here}$var_list
		rm -f ${here}_temp.txt	
	EOF
	) >> $qsub_batch
else
	# clean up
	rm -f ${here}$var_list
	rm -f ${here}_temp.txt
fi


