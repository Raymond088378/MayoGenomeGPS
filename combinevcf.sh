#!/bin/bash

if [ $# != 4 ]
then
    echo -e "script to combine multiple vcf files\nUsage: combine vcfs <input files><output vcf ><run info><to delete input files(yes/no)"
else
    set -x
    echo `date`
    input=`echo $1 | sed -e "s/ *$//" | sed -e "s/^ *//"`
    output=$2
    run_info=$3
    flag=`echo $4 | tr "[a-z]" "[A-Z]"`
    
    tool_info=$( cat $run_info | grep -w '^TOOL_INFO' | cut -d '=' -f2)
    java=$( cat $tool_info | grep -w '^JAVA' | cut -d '=' -f2)
    gatk=$( cat $tool_info | grep -w '^GATK' | cut -d '=' -f2)
    ref=$( cat $tool_info | grep -w '^REF_GENOME' | cut -d '=' -f2)
    dbSNP=$( cat $tool_info | grep -w '^dbSNP_REF' | cut -d '=' -f2)
    script_path=$( cat $tool_info | grep -w '^WORKFLOW_PATH' | cut -d '=' -f2)
	
    export PATH=$java:$PATH
    
    num_times=`echo $input | tr " " "\n" | grep -c -w "\-V"`
    if [ $num_times == 1 ]
    then
        file=`echo $input | sed -e '/-V/s///g' | sed -e "s/ *$//" | sed -e "s/^ *//"`
        cp $file $output
        if [ -s $file.idx ]
        then
            cp $file.idx $output.idx
        fi
    else    
        $java/java -Xmx2g -Xms512m -jar $gatk/GenomeAnalysisTK.jar \
        -R $ref \
        -et NO_ET \
        -K $gatk/Hossain.Asif_mayo.edu.key \
        -T CombineVariants \
        $input \
        -o $output
    fi
    if [ ! -s $output.idx ]
    then
        if [ $num_times == 1 ]
        then
            echo -e "\njust copied the files as there was only one file.\n"
            for i in `echo $input | sed -e '/-V/s///g'`
	    do
		if [ -s $i ]
		then
			rm $i
		fi
	done			
        else
            $script_path/errorlog.sh $output combinevcf.sh ERROR "failed to create"
            exit 1;
        fi
    else
        if [ $flag == "YES" ]
        then
            for i in `echo $input | tr " " "\n" | awk '$0 !~ /-V/ && $0 !~ /-priority/' | tr "\n" " "`
	    do
		if [ -s $i ]
		then
			rm $i
		fi
	done		
            sample=`echo $input | tr " " "\n" | awk '$0 !~ /-V/ && $0 !~ /-priority/' | awk '{print $0".idx"}'`
			for i in $sample
			do
				if [ $i != ".idx" ]
				then
					if [ -s $i ]
					then
						rm $i
					fi
				fi
			done	
        fi
    fi
    echo `date`
fi	
	