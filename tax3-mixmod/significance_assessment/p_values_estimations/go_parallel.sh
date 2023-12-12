#!/bin/bash

####################################################
#
# launch or sets-up job batches for selected variables and components
#
####################################################

echo -e "\n ---- mixmod setup script ----\n"

# --------------------------------------------------------------
# check parameters
#
_usage () {
cat <<- EOF
    *** go_parallel.sh ***

        Setup or launch EM algorithm jobs with gaussian or beta distribution

        <parameters_file>   list of analysis parameters sourced by this script

        --nforks=<n>    Number of bash forks used to write jobs
                        and number of qsub_batches that will be launched
                        Objective is to speed up set-up and minimize load on torque job submission
                        Default is 5

        --job-size=<n>  Number of variables assessed sequentially in each job
                        Objective is to minimize the number of qsub jobs launched
                        Default is 10
                        Consider a lower value with EM part2, ideally --job-size=1

        -h,--help       Displays this and stops

        --usage         Displays this and keeps going

EOF
}

# --------------------------------------------------------------
# get analysis options

OPTS=$(getopt -o h -l help,usage::,nforks:,job-size: -- "$@")
if [ $? != 0 ] ; then
    echo "$0 $@ getopt error"
    exit 1
fi

eval set -- "$OPTS"

# ---
# default init tax3 location and version

# number of parallel jobs (aka --nforks)
njobs=5

# size of aggregates (aka --job-size)
requestedSizeOfAggregates=10


# management of user arguments
while true
do
    case "$1" in
        -h|--help)          _usage
                            exit 0
                            ;;

        --usage)            case "$2" in
                                yes)    _usage ;;
                                no)     ;;
                                *)      _usage ;;
                            esac
                            shift # delete "$1" and "$2" after esac
                            ;;

        --nforks)           njobs="$2"
                            shift # delete "$1" and "$2" after esac
                            ;;

        --job-size)         requestedSizeOfAggregates="$2"
                            shift # delete "$1" and "$2" after esac
                            ;;

        --)                 shift # remove -- getopt tag
                            break
                            ;;

        *)                  _usage
                            echo "$0 $@  : Illegal option: $1"
                            exit 1
                            ;;
    esac
    shift  # next option
done

# parameters file is the remaining argument
if [ "$#" -eq 0 ]; then
    _usage
    die "$0 $@ : *** Parameter file missing"
elif [ "$#" -eq 1 ]; then
    USER_PARAMETER_FILE="$1"
else
    _usage
    die "$0 $@ : *** Too many argumnents: $*"
fi

# --------------------------------------------------------------
# get analysis parameters

# check whether parameters file is readable
if [ ! -r "${USER_PARAMETER_FILE}" ] ; then echo "Error: parameters file is not readable"; exit; fi

# source the parameters file
. ${USER_PARAMETER_FILE}

# --------------------------------------------------------------

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
mkdir -p $here

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

# pre run & main run modes
doPreRun=$DO_PRE_RUN
PreRunInitType=$PRE_RUN_INIT_TYPE
RunInitType=$RUN_INIT_TYPE

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

em_timeout=$EM_TIMEOUT

# step 2
enable_step2=$ENABLE_STEP2

# step2 exploration
step2_exploration_raised_limit=$STEP2_EXPLORATION_RAISED_LIMIT
step2_exploration_higher_limit=$STEP2_EXPLORATION_HIGHER_LIMIT

# are we doing graphs? (no: speed-up?)
dograph=$DOGRAPH

# title for graphs
title=$TITLE

# misc
minimize_io=$MINIMIZE_IO

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

#
# FORKED JOB
# $1 : part id
# $2 : job id offset
# $3 : data file
#
forkedjob() {

	partid=$1
	jobidoffset=$2
	datapartfile=$3

	# ----------------------------------------
	# main loop
	echo -e "Launching analyses ... $partid \n"

	# create directory where all results and working files are located
	mkdir -p $outdir

	# file listing all analyses and location (used in data management)
	case "${loading_type}" in
		component)
			echo "variable,component,jobid,path" > ${outdir}analyses${partid}.csv
			;;
		projection)
			echo "variable,dimension,jobid,path" > ${outdir}analyses${partid}.csv
			;;
	esac

	# file for qsub job removal
	rm ${here}qsub_del${partid}.sh 2> /dev/null
	touch ${here}qsub_del${partid}.sh
	chmod u+x ${here}qsub_del${partid}.sh

	# temp log file
	rm ${here}_temp${partid}.txt 2> /dev/null
	touch ${here}_temp${partid}.txt

	if [ "$qsub_batch" != "" ]
	then
		# qsub command go to a batch file
		(
		cat <<- EOF
			#!/bin/bash

			echo -e "\n Launching analyses ... \n"

			# create directory where all results and working files are located
			mkdir -p $outdir

			# local variables (aim at reducing batch size)
			here=${here}
			outdir=${outdir}

			# file listing all analyses and location (used in data management)
			case "${loading_type}" in
				component)
					echo "variable,component,jobid,path" > \${outdir}analyses${partid}.csv
					;;
				projection)
					echo "variable,dimension,jobid,path" > \${outdir}analyses${partid}.csv
					;;
			esac

			# file for qsub job removal
			rm \${here}qsub_del${partid}.sh 2> /dev/null
			touch \${here}qsub_del${partid}.sh
			chmod u+x \${here}qsub_del${partid}.sh

			# temp log file
			rm \${here}_temp${partid}.txt 2> /dev/null
			touch \${here}_temp${partid}.txt

		EOF
		) > ${here}qsub_batch${partid}.sh
		chmod u+x ${here}qsub_batch${partid}.sh
	fi

	# aggregation : jobs are pooled together in batches, to reduce number of qsub commands
	agg_nperbatch=$requestedSizeOfAggregates			# number of jobs per batch
	agg_currentnperbatch=0	# number of jobs per batch
	agg_n=0						# number of aggregated jobs

	# for each (var,comp)
	ijob=${jobidoffset}
	nlines=`wc -l ${datapartfile} | cut --delimiter=' ' --field=1`
	while read varcomp ;
	do
		ijob=$((1+$ijob))
		nlines=$(($nlines - 1))

		# get varname and component to analyse
		var=`echo $varcomp | cut --delimiter=, --fields=1 | sed -e"s/ //g"`
		comp=`echo $varcomp | cut --delimiter=, --fields=2 | sed -e"s/ //g"`

		# where to store results
		case "${loading_type}" in
			component)
				vardir="/comp"$comp"_"$var"/"
				# log
				echo -ne "${ijob}  "
				printf "%15s" $var
				budimcomp='", component " $comp "   ...   "'
				echo -ne ", component " $comp "   ...   \n"
				;;
			projection)
				vardir="/dim"$comp"_"$var"/"
				# log
				echo -ne "${ijob}  "
				printf "%15s" $var
				budimcomp='", dimension " $comp "   ...   "'
				echo -ne ", dimension " $comp "   ...   \n"
				;;
		esac

		#
		# make working directory
		# we have to use $partid here to reduce the number of directories in $outdir (limit is 32K)
		resdir=${outdir}${partid}${vardir}results
		workdir=${outdir}${partid}${vardir}workdir

		# first remove existing directory
		rm -fr $workdir $resdir
		# then create working directory
		mkdir -p $workdir $resdir

		#
		# copy header and observed loading
		#
		if [ -f $loadings ] ; then
			# a file with all observed loadings : have to grep, and get proper column
			# csv file, spaces allowed
			# "rs1234 ,0.5"
			head -1 $loadings    | cut --delimiter=, --fields=1,$((1+$comp)) >  ${resdir}/loading.csv
			grep "^${var}\( \)\?\," $loadings  | cut --delimiter=, --fields=1,$((1+$comp)) >> ${resdir}/loading.csv
		else
			# a directory with one file for each variable (much faster) : just have to get proper column
			cat ${loadings}/${var}.csv  | cut --delimiter=, --fields=1,$((1+$comp)) >  ${resdir}/loading.csv
		fi

		# check we have only 2 lines:
		skip=false
		if [[ `wc -l ${resdir}/loading.csv | cut --fields=1 --delimiter=' '` != '2' ]] ; then
		   echo "error: loading file does not have 2 lines, but `wc -l ${resdir}/loading.csv | cut --fields=1 --delimiter=' '` : "
		   cat ${resdir}/loading.csv
		   echo "=> skipping"
		   skip=true
		fi

		if [ $skip = false ] ; then
			# copy resamplings
			n=$((1+$n_resamplings))
			compcol=$(( $comp + 2 ))  # no p & CI in C++
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
					-e"s/<DO_PRE_RUN>/"$doPreRun"/"             \
					-e"s/<PRE_RUN_INIT_TYPE>/"$PreRunInitType"/"             \
					-e"s/<RUN_INIT_TYPE>/"$RunInitType"/"             \
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
					-e"s/<DOGRAPH>/$dograph/"                 \
					-e"s/<MINIMIZE_IO>/${minimize_io}/"                 \
					-e"s/<RUN_CHAINS_IN_PARALLEL>/${RunChainsInParallel}/"                 \
					-e"s/<TITLE>/$title/"                 \
					-e"s/<PCA_COMPONENT>/"$comp"/"        \
					\
					${batchdir}$runner > $workdir/$runner


			chmod u+x $workdir/$runner

			# EM algorithm needed if beta model
			cp ${batchdir}EM_functions.R $workdir; chmod u+x $workdir/EM_functions.R

			# launch batch job
			if [ "$exec_mode" = "qsub" ]
			then
				if [ "${qsub_batch}${partid}" = "" ]
				then
					# qsub
					qsub $qsub_args -N ${jobname} -o $workdir -e $workdir $workdir/$runner | tee ${here}_temp${partid}.txt
					echo -n "qdel " >> ${here}qsub_del${partid}.sh
					cat ${here}_temp${partid}.txt >> ${here}qsub_del${partid}.sh
					echo >> ${here}qsub_del${partid}.sh
					# log
					pbsjobid=`head -1 ${here}_temp${partid}.txt`
					echo ${var},${comp},${pbsjobid},${resdir} >> ${outdir}analyses${partid}.csv
				else
					if [ ${agg_currentnperbatch} == 0 ] ; then
						# create new aggregate
						agg_n=$(( $agg_n + 1 ))
						(
						cat <<- EOF
							#!/bin/bash

							outdir=${outdir}
							aggdir=\${outdir}${partid}

						EOF
						) 	> ${outdir}${partid}/part${partid}_agg${agg_n}.sh

					fi

					# add new entry to aggregate
					agg_currentnperbatch=$(( ${agg_currentnperbatch} + 1 ))
					(
					cat <<- EOF

						cd	\${aggdir}${vardir}workdir
						timeout ${em_timeout} ./${runner}
						echo ${var},${comp},\${PBS_JOBID},${resdir} >> ${outdir}analyses${partid}.csv

					EOF
					) 	>> ${outdir}${partid}/part${partid}_agg${agg_n}.sh

					# if aggregate full or last line , qsub it
					if [  ${agg_currentnperbatch} -ge ${agg_nperbatch} ] || [ $nlines == 0 ] ; then
						# reset
						agg_currentnperbatch=0
						# qsub command go to a batch file
						(
						cat <<- EOF
							# submit new job
							workdir=${outdir}${partid}
							resdir=${outdir}${partid}
							echo -ne "${agg_n}  "
							qsub $qsub_args -N ${jobname} -o \${workdir} -e \${workdir} \${workdir}/part${partid}_agg${agg_n}.sh | tee \${here}_temp${partid}.txt
							echo -n "qdel " >> \${here}qsub_del${partid}.sh
							cat \${here}_temp${partid}.txt >> \${here}qsub_del${partid}.sh
							echo >> \${here}qsub_del${partid}.sh
							sleep 0.5

						EOF
						) >> ${here}qsub_batch${partid}.sh
					fi
				fi
			else
				# at command
				at -f $workdir/$runner now  | tee ${here}_temp${partid}.txt
				# log
				jobid=`head -1 ${here}_temp${partid}.txt`
				echo $var","$comp","$jobid","$resdir >> ${outdir}analyses.csv
			fi

		fi # skipping

	done < ${datapartfile} # next variable

	# clean up
	if [ "$exec_mode" = "qsub" ] && [ "$qsub_batch" != "" ]
	then
		(
		cat <<- EOF
			# clean up
			rm -f ${datapartfile}
			rm -f ${here}_temp${partid}.txt
		EOF
		) >> ${here}qsub_batch${partid}.sh
	else
		# clean up
		rm -f ${datapartfile}
		rm -f ${here}_temp${partid}.txt
	fi

} # end of forkedjob function


# ----------------------------------------
# cut the header from list of variables
# make sure there is no (var,comp) duplicate
cat $variables | sed '1d' | sort | uniq > ${here}$var_list

# define and store jobname for job retrieval
jobname="EM_${RANDOM}"
echo ${jobname} > ${here}/jobname.txt

#
# split list of variables and
# run multiple instances of forkedjob main function on each part
#
rm ${here}${var_list}_part.* 2> /dev/null

# get number of lines ans split accordingly to njobs
nlines=$( wc -l ${here}$var_list | cut --field=1 --delimiter=" " )
nlinespart=$(( 1 + $nlines / $njobs ))
split -l $nlinespart ${here}$var_list ${here}${var_list}_part.

nparts=0
jobidstart=0
allparts=$( ls ${here}${var_list}_part.* )
for part in $allparts ; do
    nparts=$(( $nparts + 1 ))

    # launch this part
    forkedjob ${nparts} ${jobidstart} $part &

    # get its pid
    pids[$nparts]=$!

    # increase jobid offset
    jobidstart=$(( ${jobidstart} + ${nlinespart} ))

done

# wait for all jobs to end
echo "waiting for completion of $nparts jobs:"
for i in $( seq 1 $nparts ) ; do
    wait ${pids[$i]}
done

#
# if qsub batch : create batch calling all sub-batches
#
if [ "$exec_mode" = "qsub" ] && [ "$qsub_batch" != "" ]
then
    cat > $qsub_batch <<- EOF
		#!/bin/bash

		#
		# batch created automagically by go_parallel.sh script
		#

		# local variables
		here=${here}
		outdir=${outdir}

		# call sub-batches
		for i in \$( seq 1 $nparts ) ; do
		    \${here}qsub_batch\${i}.sh
		done

		# merge analyses*.csv once all jobs submitted
		head -1 \${outdir}analyses1.csv > \${outdir}analyses.csv
		for i in \$( seq 1 $nparts ) ; do
		    tail -n+2 \${outdir}analyses\${i}.csv >> \${outdir}analyses.csv
		done
	EOF
    chmod u+x $qsub_batch

else
    head -1 ${outdir}analyses1.csv > ${outdir}analyses.csv
    for i in $( seq 1 $nparts ) ; do
        tail -n+2 ${outdir}analyses${i}.csv >> ${outdir}analyses.csv
    done
fi

#
# remove all jobs batch
#
cat > ./qdel_all.sh <<- EOF
	#!/bin/bash

	#
	# batch created automagically by go_parallel.sh script
	#

	# local variables
	here=${here}

	# kill all sub-batches (in reverse order)
	for i in \$( seq 1 ${nparts} | tac ) ; do
	    sort -r \${here}qsub_del\${i}.sh > _temp_del.sh
	    chmod u+x _temp_del.sh
	    ./_temp_del.sh
	done

	# WARNING !!! kill all EM processes and Torque jobs in all user hosts
	if true ; then
	    cat /etc/hosts | while read host
	    do
	        if [[ "\${host}" == *"\${USER}"* ]]; then
	            dns=\$( echo \${host} | cut --delimiter=' ' --fields=3 )
	            echo \${dns}

	            RET=\$(
	                ssh -T \${dns} <<- ENDSSH
	                    pkill -f .fileserver.
					ENDSSH
	                )

	            RET=\$(
	                ssh -T \${dns} <<- ENDSSH
	                    pkill -f runner_mixture.sh
					ENDSSH
	                )

	        fi
	    done
	fi
EOF

chmod u+x ./qdel_all.sh
