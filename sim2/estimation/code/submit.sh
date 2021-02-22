#!/bin/bash

# when you execute this script on the cluster run this command
# chmod u+x submit.sh
# ./submit.sh run_sim.R simulation

##################### Change these constants ##############################
analysis=$2      # change for every analysis you run (2nd arg)
maildom='@emory.edu'   # your email domain (for receiving error messages)
myscratch="path/scratch"  # location of your persistent scratch dir
resultdir="path/out"  # This is a folder in permanent storage
script=$1      # your code as (R or Python) script (1st arg)
max_jobs=10   # max number of jobs to run at a time
total_jobs=30 # needs to be divisible by max_jobs
############## typically you don't have to change anything below here #######

username=$(id -nu)

# run R script to get list size
# this apparently doesn't work
# listsize=$(${script} listsize)

nloops=$((${total_jobs}/${max_jobs}-1))

# if scratch directory doesn't exist, make it
[ ! -d ${myscratch} ] && mkdir ${myscratch}
[ ! -d ${myscratch}/out ] && mkdir ${myscratch}/out
[ ! -d ${myscratch}/err ] && mkdir ${myscratch}/err

# submit first batch of jobs
for i in $(seq 1 ${max_jobs}); do
	echo "#!/bin/bash" >> script$i.out
	echo "#SBATCH --nodes=1 # ask for 1 node" >> script$i.out
	echo "#SBATCH --ntasks-per-node=10 # 4 tasks each node" >> script$i.out
	echo "#SBATCH --mem-per-cpu=1G" >> script$i.out
	echo "#SBATCH --partition=guo-cbis" >> script$i.out
	echo "#SBATCH --job-name=${analysis}$i" >> script$i.out
	echo "#SBATCH --error=${myscratch}/err/${analysis}$i.err" >> script$i.out
	echo "#SBATCH --output=${myscratch}/out/${analysis}$i.out" >> script$i.out

	echo "module load R/4.0.2" >> script$i.out
	echo "R CMD BATCH '--args iter=$i' ${script}" >> script$i.out

	chmod +x script$i.out
	sbatch script$i.out
done

# sub next batches of jobs holding 
for j in $(seq 1 $nloops); do
	for i in $(seq 1 ${max_jobs}); do
		hold_name=${analysis}$((($j-1)*${max_jobs}+$i))
		jid=$(($j*${max_jobs}+$i))

		echo "#!/bin/bash" >> script$jid.out
		echo "#SBATCH --nodes=1 # ask for 1 node" >> script$jid.out
		echo "#SBATCH --ntasks-per-node=10 # 4 tasks each node" >> script$jid.out
		echo "#SBATCH --mem-per-cpu=1G" >> script$jid.out
		echo "#SBATCH --partition=guo-cbis" >> script$jid.out
		echo "#SBATCH --job-name=${analysis}$jid" >> script$jid.out
		echo "#SBATCH --error=${myscratch}/err/${analysis}$jid.err" >> script$jid.out
		echo "#SBATCH --output=${myscratch}/out/${analysis}$jid.out" >> script$jid.out


		echo "module load R/4.0.2" >> script$jid.out
		echo "R CMD BATCH '--args iter=$jid' ${script}" >> script$jid.out

		chmod +x script$jid.out
		sbatch script$jid.out
	done
done

