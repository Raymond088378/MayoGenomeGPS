#!/bin/bash
if [ $# != 1 ]
then
	echo "Usage <limit of jobs>"
else
	limit=$1
	sleep 10
	let count=1
	num_jobs=`qstat | awk '{print $1}' | grep -E '^[0-9]+$' | sort | uniq | wc -l`
	while [ $num_jobs -ge $limit ]
	do
		let wait=`expr 600 "*" $count`
		echo "waiting for slot on cluster"
		sleep $wait
		num_jobs=`qstat | awk '{print $1}' | grep -E '^[0-9]+$' | sort | uniq | wc -l`
		let count=count+1
	done
fi	