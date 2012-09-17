#!/bin/bash
if [ $# != 1 ]
then
	echo -e "script to check the number of jobs a user can submit\nUsage: <limit of jobs>"
else
	limit=$1
	sleep 10
	let count=1
	TO=`id |awk -F '(' '{print $2}' | cut -f1 -d ')'`
	num_jobs=`qstat -u $TO | awk '{print $1}' | grep -E '^[0-9]+$' | sort | uniq | wc -l`
	while [ $num_jobs -ge $limit ]
	do
		if [ $count -eq 1 ]
		then
			echo -e "User reached the limit of $limit jobs on RCF cluster so the submission script is waiting for free slots" |  mailx -v -s "JOB LIMIT on RCF cluster" "$TO" 
		fi	
		let wait=`expr 600 "*" $count`
		echo "waiting for slot on cluster"
		sleep $wait
		num_jobs=`qstat | awk '{print $1}' | grep -E '^[0-9]+$' | sort | uniq | wc -l`
		let count=count+1
	done
	jobs=`expr $limit "-" $num_jobs`
	echo "You can submit $jobs more jobs" 
fi	