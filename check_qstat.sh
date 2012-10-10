#!/bin/bash

if [ $# != 1 ]
then
	echo -e "script to check the number of jobs a user can submit\nUsage: ./check_qstat.sh <limit of jobs>"
else
    limit=$1
	sleep 10
	let count=1
	num=`qstat | awk '{ if ($1>0) print $1;}'` 
	if [[ "$num" ]]
	then
		num_jobs=`echo $num | tr " " "\n" | wc -l `
	else
		num_jobs=$limit
	fi	
	while [ $num_jobs -ge $limit ]
	do
		if [ $count -eq 1 ]
		then
			echo -e "User reached the limit of $limit jobs on RCF cluster so the submission script is waiting for free slots" |  mailx -v -s "JOB LIMIT on RCF cluster" "$TO" 
		fi	
		let wait=`expr 60 "*" $count`
		echo "waiting for slot on cluster"
		sleep $wait
		num=`qstat | awk '{ if ($1>0) print $1;}'`
		if [[ "$num" ]]
		then
			num_jobs=`echo $num | tr " " "\n" | wc -l `
		else
			num_jobs=$limit
		fi
		let count=count+1
	done
	jobs=`expr $limit "-" $num_jobs`
	echo "You can submit $jobs more jobs" 
fi	