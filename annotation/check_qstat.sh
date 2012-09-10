#!/bin/bash

num_jobs=`qstat | awk '{print $1}' | sort | uniq | wc -l`

while [ $num_jobs -ge 2500 ]
do
	sleep 5m
	num_jobs=`qstat | awk '{print $1}' | sort | uniq | wc -l`
done
	