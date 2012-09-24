#!/bin/bash

if [ $# != 4 ]
then
    echo -e " Usage: script to merge the covergae pileup \n <input dir><output dir><sample><run info >"
else
    set -x
    echo `date`
    input=$1
    output=$2
    sample=$3
	run_info=$4
   
    ### expected file name for covergae pileup file is ${sample}.chr${chr}.pileup.i.out
    tool_info=$( cat $run_info | grep -w '^TOOL_INFO' | cut -d '=' -f2 | tr ":" " " )
    chrs=$( cat $run_info | grep -w '^CHRINDEX' | cut -d '=' -f2 | tr ":" " " )
    multi=$( cat $run_info | grep -w '^MULTISAMPLE' | cut -d '=' -f2| tr "[a-z]" "[A-Z]")
    sample_info=$( cat $run_info | grep -w '^SAMPLE_INFO' | cut -d '=' -f2)
	script_path=$( cat $tool_info | grep -w '^WORKFLOW_PATH' | cut -d '=' -f2 )
    #### assuming that each covergae file has 40 rows
    
    if [ $multi == "YES" ]
    then
        pair=$( cat $sample_info | grep -w "^$sample" | cut -d '=' -f2)
        for i in $pair
        do
            for chr in $chrs
            do
                num_rows=`cat $input/$sample.$i.chr$chr.pileup.i.out | wc -l`
                if [ $num_rows != 100 ]
                then
                    $script_path/errorlog.sh $input/$sample.$i.chr$chr.pileup.i.out getCoverage.sh ERROR "malformed"
					exit 1;
                fi
            done	

            for ((k=1;k<=100;k++));
            do
                total=0
                for j in $chrs
                do
                    a=`sed -n "$k{p;q}" $input/$sample.$i.chr$j.pileup.i.out`
                    total=`expr $total "+" $a`
                done
                echo $total >> $output/$sample.$i.coverage.out
            done
        done
    else
        for chr in $chrs
        do
            num_rows=`cat $input/$sample.chr$chr.pileup.i.out | wc -l`
            if [ $num_rows != 100 ]
            then
                $script_path/errorlog.sh $input/$sample.chr$chr.pileup.i.out getCoverage.sh ERROR "malformed"
				exit 1;
            fi
        done	
        
        for ((i=1;i<=100;i++));
        do
            total=0
            for j in $chrs
            do
                a=`sed -n "$i{p;q}" $input/$sample.chr$j.pileup.i.out`
                total=`expr $total "+" $a`
            done
            echo $total >> $output/$sample.coverage.out
        done
    fi
    echo `date`
fi