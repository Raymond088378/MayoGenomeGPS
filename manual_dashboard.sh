#!/bin/bash

if [ $# != 4 ]
then
	echo -e "script to upload the status to the dash board when megng diffrent flowcell\nUsage: ./manual_dashboard.sh <','seperated flowcells> <status <Beginning/Complete/Delivered> <type of sample (Exome/WholeGenome)><full/path/to/run info file>"
	exit 1;
fi	
flowcells=`echo $1 | tr "," " "`
status=$2
tool=$3
run_info=$4
tool_info=$( cat $run_info | grep -w '^TOOL_INFO' | cut -d '=' -f2)
script_path=$( cat $tool_info | grep -w '^WORKFLOW_PATH' | cut -d '=' -f2)
samples=$( cat $run_info | grep -w '^SAMPLENAMES' | cut -d '=' -f2)

for flow in $flowcells
do
	id=`echo $flow |awk -F'_' '{print $NF}' | sed 's/.\(.*\)/\1/'`
	for sample in `echo $samples | tr ":" " "`
	do
		pos=$( cat $run_info | grep -w '^SAMPLENAMES' | cut -d '=' -f2 | tr ":" "\n" | grep -w -n $sam | cut -d ":" -f1)
		lanes=$( cat $run_info | grep -w '^LANEINDEX' | cut -d '=' -f2 | tr ":" "\n" | head -n $pos | tail -n 1 | tr "," " ")
        i=1
        for lane in $lanes
        do
            index=$( cat $run_info | grep -w '^LABINDEXES' | cut -d '=' -f2 | tr ":" "\n" | head -n $pos | tail -n 1 | tr "," "\n" | head -n $i | tail -n 1)
			if [ $index == "-" ]
			then
				$java/java $mem -jar $script_path/AddSecondaryAnalysis.jar -p $script_path/AddSecondaryAnalysis.properties -c -l $lane -f $id -r $flow -s $status -a $tool
			else
				$java/java $mem -jar $script_path/AddSecondaryAnalysis.jar -p $script_path/AddSecondaryAnalysis.properties -c -l $lane -i $index -f $id -r $flow -s $status -a $tool
			fi
		done
	done
done
	