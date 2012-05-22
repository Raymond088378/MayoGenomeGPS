#!/bin/sh

if [ $# != 6 ]
then
	echo -e "Usage: script to run snvmix \n <sample> <pileup file><output vcf> <mode><chr><run info>"
else
    set -x
    echo `date`
    
    sample=$1
    pileup=$2
    output=$3
    mode=`echo $4 | tr "[A-Z]" "[a-z]"`
    kit=$5
    run_info=$6
    
    tool_info=$( cat $run_info | grep -w '^TOOL_INFO' | cut -d '=' -f2)
    snvmix=$( cat $tool_info | grep -w '^SNVmix' | cut -d '=' -f2)
    only_ontarget=$( cat $tool_info | grep -w '^TARGETTED' | cut -d '=' -f2 | tr "[a-z]" "[A-Z]" )
    tool=$( cat $run_info | grep -w '^TYPE' | cut -d '=' -f2|tr "[A-Z]" "[a-z]")
    script_path=$( cat $tool_info | grep -w '^WHOLEGENOME_PATH' | cut -d '=' -f2 )
    prob_filter=$( cat $tool_info | grep -w '^PROB_FILTER' | cut -d '=' -f2)
    TargetKit=$( cat $tool_info | grep -w '^ONTARGET' | cut -d '=' -f2 )
    bedtools=$( cat $tool_info | grep -w '^BEDTOOLS' | cut -d '=' -f2 )
            
    temp=`echo $output | sed -e '/.vcf/s///g'`
    temp_kit=`echo $kit | sed -e '/-L/s///g' | sed -e "s/^ *//"`
    if [ $mode == "all" ]
    then	
        $snvmix/SNVMix2 -i $pileup -f -m $snvmix/Mu_pi.txt -o $temp	
    else
        $snvmix/SNVMix2 -i $pileup -m $snvmix/Mu_pi.txt -o $temp
    fi

    if [[ $only_ontarget == "YES" && $tool == "exome" ]]
    then
        perl $script_path/snvmix_to_vcf.pl -i $temp -o $output -s $sample -p $prob_filter
        if [ -s $temp_kit ]
		then
            $bedtools/intersectBed -a $output -b $temp_kit -wa -header > $output.i
			rm $temp_kit
        else
            cp $output $output.i	
        fi	
        mv $output.i $output
    else
        perl $script_path/snvmix_to_vcf.pl -i $temp -o $output -s $sample
    fi
    
    if [ -s $output ]
    then
        rm $temp
        rm $pileup
    else
		$script_path/errorlog.sh $output snvmix2.sh ERROR "failed to create"	
        exit 1;
    fi	
    echo `date`
fi	
		