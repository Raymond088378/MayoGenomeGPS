#!/bin/bash

if [ $# -le 3 ]
then
    echo "Usage : script to update secondary dashboard \n <sample ><runinfo ><stage of the workflow> <status of the stage> <id>"
else
    echo `date`	
    sample=$1
    run_info=$2
    stage=$3
    status=$4
    if [ $5 ]
    then
        id=$5
    fi
    analysis=$( cat $run_info | grep -w '^ANALYSIS' | cut -d '=' -f2 |tr "[A-Z]" "[a-z]")
    tool_info=$( cat $run_info | grep -w '^TOOL_INFO' | cut -d '=' -f2)
    script_path=$( cat $tool_info | grep -w '^WORKFLOW_PATH' | cut -d '=' -f2 )
    java=$( cat $tool_info | grep -w '^JAVA' | cut -d '=' -f2)
    run_num=$( cat $run_info | grep -w '^OUTPUT_FOLDER' | cut -d '=' -f2)
    flowcell=`echo $run_num | awk -F'_' '{print $NF}' | sed 's/.\(.*\)/\1/'`
    version=$( cat $run_info | grep -w '^VERSION' | cut -d '=' -f2)
	
    
	##tool is hard coded
	type=$( cat $run_info | grep -w '^TYPE' | cut -d '=' -f2| tr "[A-Z]" "[a-z]")
	if [ $type == "exome" ]
	then
		tool=Exome
	else
		tool=WholeGenome
    fi
	
    if [ $analysis == "mayo" -o $analysis == "realign-mayo" ]
    then
        if [ $5 ]
        then
            pos=$( cat $run_info | grep -w '^SAMPLENAMES' | cut -d '=' -f2 | tr ":" "\n" | grep -n $sample | cut -d ":" -f1)
            lane=$( cat $run_info | grep -w '^LANEINDEX' | cut -d '=' -f2 | tr ":" "\n" | head -n $pos | tr "," "\n" | head -n $id | tail -n 1)
            index=$( cat $run_info | grep -w '^LABINDEXES' | cut -d '=' -f2 | tr ":" "\n" | head -n $pos | tr "," "\n" | head -n $id | tail -n 1) 
            if [ $index == "-" ]
            then
                if [ $status == "complete" ]
                then
                    $java/java -Xmx32M -jar $script_path/AddSecondaryAnalysis.jar -p $script_path/AddSecondaryAnalysis.properties -c -l $lane -f $flowcell -r $run_num -s $stage -a $tool
                else
                    $java/java -Xmx32M -jar $script_path/AddSecondaryAnalysis.jar -p $script_path/AddSecondaryAnalysis.properties -l $lane -f $flowcell -r $run_num -s $stage -a $tool
                fi
            else
                if [ $status == "complete" ]
                then
                    $java/java -Xmx32M -jar $script_path/AddSecondaryAnalysis.jar -p $script_path/AddSecondaryAnalysis.properties -c -l $lane -f $flowcell -i $index -r $run_num -s $stage -a $tool
                else
                    $java/java -Xmx32M -jar $script_path/AddSecondaryAnalysis.jar -p $script_path/AddSecondaryAnalysis.properties -l $lane -f $flowcell -i $index -r $run_num -s $stage -a $tool
                fi		
            fi
        else
            s=`echo $sample | tr ":" " "`
            for sam in $s
            do
                pos=$( cat $run_info | grep -w '^SAMPLENAMES' | cut -d '=' -f2 | tr ":" "\n" | grep -n $sam | cut -d ":" -f1)
                lanes=$( cat $run_info | grep -w '^LANEINDEX' | cut -d '=' -f2 | tr ":" "\n" | head -n $pos | tail -n 1 | tr "," " ")
                i=1
                for lane in $lanes
                do
                    index=$( cat $run_info | grep -w '^LABINDEXES' | cut -d '=' -f2 | tr ":" "\n" | head -n $pos | tail -n 1 | tr "," "\n" | head -n $i | tail -n 1)
                    if [ $index == "-" ]
                    then
                        if [ $status == "complete" ]
                        then
                            $java/java -Xmx32M -jar $script_path/AddSecondaryAnalysis.jar -p $script_path/AddSecondaryAnalysis.properties -c -l $lane -f $flowcell -r $run_num -s $stage -a $tool
                        else
                            $java/java -Xmx32M -jar $script_path/AddSecondaryAnalysis.jar -p $script_path/AddSecondaryAnalysis.properties -l $lane -f $flowcell -r $run_num -s $stage -a $tool
                        fi		
                    else
                        if [ $status == "complete" ]
                        then
                            $java/java -Xmx32M -jar $script_path/AddSecondaryAnalysis.jar -p $script_path/AddSecondaryAnalysis.properties -c -l $lane -f $flowcell -i $index -r $run_num -s $stage -a $tool
                        else
                            $java/java -Xmx32M -jar $script_path/AddSecondaryAnalysis.jar -p $script_path/AddSecondaryAnalysis.properties -l $lane -f $flowcell -i $index -r $run_num -s $stage -a $tool
                        fi		
                    fi
                    let i=i+1
                done		
            done
        fi
    fi    
    echo `date`
fi


