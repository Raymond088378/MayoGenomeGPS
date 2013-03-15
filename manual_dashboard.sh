#!/bin/bash

if [ $# != 5 ]
then
	echo -e "script to upload the status to the dash board when merging different flow cell\nUsage: ./manual_dashboard.sh <flowcell/run ID> <status <Beginning/Complete/Delivered> <type of sample (Exome/WholeGenome)><full/path/to/run info file> <date of the analysis (mm/dd/yy)>"
	exit 1;
fi
flowcells=$1
status=$2
tool=$3
run_info=$4
date=$5
tool_info=$( cat $run_info | grep -w '^TOOL_INFO' | cut -d '=' -f2)
script_path=$( cat $tool_info | grep -w '^WORKFLOW_PATH' | cut -d '=' -f2)
samples=$( cat $run_info | grep -w '^SAMPLENAMES' | cut -d '=' -f2)
java=$( cat $tool_info | grep -w '^JAVA' | cut -d '=' -f2)
for flow in $flowcells
do
	id=`echo $flow |awk -F'_' '{print $NF}' | sed 's/.\(.*\)/\1/'`
	for sample in `echo $samples | tr ":" " "`
	do
		pos=$( cat $run_info | grep -w '^SAMPLENAMES' | cut -d '=' -f2 | tr ":" "\n" | grep -w -n $sample | cut -d ":" -f1)
		lanes=$( cat $run_info | grep -w '^LANEINDEX' | cut -d '=' -f2 | tr ":" "\n" | head -n $pos | tail -n 1 | tr "," " ")
        i=1
        for lane in $lanes
        do
            index=$( cat $run_info | grep -w '^LABINDEXES' | cut -d '=' -f2 | tr ":" "\n" | head -n $pos | tail -n 1 | tr "," "\n" | head -n $i | tail -n 1)
			if [ $index == "-" ]
			then
				$java/java $mem -jar $script_path/AddSecondaryAnalysis.jar -p $script_path/AddSecondaryAnalysis.properties -c -l $lane -f $id -r $flow -s $status -a $tool -d $date
			else
				$java/java $mem -jar $script_path/AddSecondaryAnalysis.jar -p $script_path/AddSecondaryAnalysis.properties -c -l $lane -i $index -f $id -r $flow -s $status -d $date -a $tool
			fi
		done
	done
done
	