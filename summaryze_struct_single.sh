#!/bin/bash

if [ $# -le 2 ]
then	
	echo -e "script to merge the structural varaints from multiple tools\nUsage: <sample name> <base dir> <path to run_info file> <group optional> ";
	exit 1
else
	set -x
	echo `date`
	sam=$1
	basedir=$2
	run_info=$3
	if [ $4 ]
	then
		sample=$4
	else
		sample=$sam
	fi	
	
	chrs=$( cat $run_info | grep -w '^CHRINDEX' | cut -d '=' -f2 | tr ":" "\n" )
	tool_info=$( cat $run_info | grep -w '^TOOL_INFO' | cut -d '=' -f2 | tr ":" "\n" )
	output=$( cat $run_info | grep -w '^BASE_OUTPUT_DIR' | cut -d '=' -f2)
	PI=$( cat $run_info | grep -w '^PI' | cut -d '=' -f2)
	java=$( cat $tool_info | grep -w '^JAVA' | cut -d '=' -f2)
	gatk=$( cat $tool_info | grep -w '^GATK' | cut -d '=' -f2)
	tool=$( cat $run_info | grep -w '^TYPE' | cut -d '=' -f2 | tr "[A-Z]" "[a-z]" )
	run_num=$( cat $run_info | grep -w '^OUTPUT_FOLDER' | cut -d '=' -f2)
	ref=$( cat $tool_info | grep -w '^REF_GENOME' | cut -d '=' -f2 )
	script_path=$( cat $tool_info | grep -w '^WORKFLOW_PATH' | cut -d '=' -f2 )

	mkdir -p $basedir/Reports_per_Sample
	output=$basedir/Reports_per_Sample
	mkdir -p $output/SV

	#Summaryzing CNVs


	inputargs=""
	inputargs_filter=""
	input=""
	for chr in $chrs
	do
		inputfile=$basedir/cnv/$sample/$sam.$chr.del.vcf
		input=$basedir/cnv/$sample/$sam.$chr.del.bed
		if [ ! -f $inputfile ]
		then
			touch $inputfile.fix.log
			$script_path/email.sh $inputfile "not exist" run_cnvnator.sh $run_info
			$script_path/wait.sh $inputfile.fix.log 
			inputargs="-V $inputfile "$inputargs  
			cat $input >> $basedir/cnv/$sample/$sam.cnv.bed
			rm $input
		else
			inputargs="-V $inputfile "$inputargs  
			cat $input >> $basedir/cnv/$sample/$sam.cnv.bed
			rm $input
		fi
		
		inputfile=$basedir/cnv/$sample/$sam.$chr.dup.vcf 
		input=$basedir/cnv/$sample/$sam.$chr.dup.bed
		if [ ! -f $inputfile ]
		then	
			touch $inputfile.fix.log
			$script_path/email.sh $inputfile "not exist" run_cnvnator.sh $run_info
			$script_path/wait.sh $inputfile.fix.log 
			inputargs="-V $inputfile "$inputargs  
			cat $input >> $basedir/cnv/$sample/$sam.cnv.bed
			rm $input
		else
			inputargs="-V $inputfile "$inputargs  
			cat $input >> $basedir/cnv/$sample/$sam.cnv.bed
			rm $input
		fi
		
		inputfile=$basedir/cnv/$sample/$sam.$chr.final.del.vcf
		input=$basedir/cnv/$sample/$sam.$chr.final.del.bed
		if [ ! -f $inputfile ]
		then	
			touch $inputfile.fix.log
			$script_path/email.sh $inputfile "not exist" run_cnvnator.sh $run_info
			$script_path/wait.sh $inputfile.fix.log 
			cat $input >> $basedir/cnv/$sample/$sam.cnv.final.bed
			rm $input
			inputargs_filter="-V $inputfile "$inputargs_filter  
		else
			cat $input >> $basedir/cnv/$sample/$sam.cnv.final.bed
			rm $input
			inputargs_filter="-V $inputfile "$inputargs_filter  
		fi
		
		inputfile=$basedir/cnv/$sample/$sam.$chr.final.dup.vcf
		input=$basedir/cnv/$sample/$sam.$chr.final.dup.bed
		if [ ! -f $inputfile ]
		then
			touch $inputfile.fix.log
			$script_path/email.sh $inputfile "not exist" run_cnvnator.sh $run_info
			$script_path/wait.sh $inputfile.fix.log 
			cat $input >> $basedir/cnv/$sample/$sam.cnv.final.bed
			rm $input
			inputargs_filter="-V $inputfile "$inputargs_filter 
		else
			cat $input >> $basedir/cnv/$sample/$sam.cnv.final.bed
			rm $input
			inputargs_filter="-V $inputfile "$inputargs_filter 
		fi
	done

	$script_path/combinevcf.sh "$inputargs" $output/SV/$sam.cnv.vcf $run_info yes
	$script_path/combinevcf.sh "$inputargs_filter" $output/SV/$sam.cnv.final.vcf $run_info yes
	
	
   ## Summaryzing Breakdancer
	inputargs=""
	input=""
	for chr in $chrs
	do
		inputfile=$basedir/struct/break/$sample/$sam.$chr.break.vcf 
		input=$basedir/struct/break/$sample/$sam.$chr.break
		if [ ! -s $inputfile ]
		then      
			touch $inputfile.fix.log
			$script_path/email.sh $inputfile "not exist" run_breakdancer.sh $run_info
			$script_path/wait.sh $inputfile.fix.log 
			cat $inputfile | awk '$0 ~/^#/' > $basedir/struct/break/$sample/$sam.header.break
			cat $inputfile | awk '$0 !~ /^#/' >> $output/SV/$sam.break.vcf
			cat $input | awk '$0 !~ /^#/' >> $basedir/struct/break/$sample/$sam.break
			rm $input
			rm $inputfile
		else
			inputargs="-V $inputfile "$inputargs
			cat $inputfile | awk '$0 ~/^#/' > $basedir/struct/break/$sample/$sam.header.break
			cat $inputfile | awk '$0 !~ /^#/' >> $output/SV/$sam.break.vcf
			cat $input | awk '$0 !~ /^#/' >> $basedir/struct/break/$sample/$sam.break
			rm $input
			rm $inputfile
		fi
	done
    cat $basedir/struct/break/$sample/$sam.inter.break | awk '$0 !~ /^#/' >> $basedir/struct/break/$sample/$sam.break
    cat $basedir/struct/break/$sample/$sam.inter.break.vcf |  awk '$0 !~ /^#/' >> $output/SV/$sam.break.vcf
	rm 	$basedir/struct/break/$sample/$sam.inter.break.vcf
	rm $basedir/struct/break/$sample/$sam.inter.break 
    cat $basedir/struct/break/$sample/$sam.header.break $output/SV/$sam.break.vcf > $output/SV/$sam.break.vcf.temp
	mv $output/SV/$sam.break.vcf.temp $output/SV/$sam.break.vcf
	rm $basedir/struct/break/$sample/$sam.header.break

	if [ ! -s $output/SV/$sam.break.vcf ]
	then
		$script_path/errorlog.sh $output/$sam.break.vcf summaryze_struct_single.sh ERROR "not created"
		exit 1;
	fi
	
	#Summaryzing Crest
	inputargs=""
    inputargs_filter=""
    input=""
    input_filter=""
    for chr in $chrs
    do
        inputfile=$basedir/struct/crest/$sample/$sam.$chr.vcf
        inputfile_filter=$basedir/struct/crest/$sample/$sam.$chr.final.vcf
        input=$basedir/struct/crest/$sample/$sam.$chr.predSV.txt
        input_filter=$basedir/struct/crest/$sample/$sam.$chr.final.predSV.txt
        if [ ! -s $inputfile ]
        then      
            touch $inputfile.fix.log
			$script_path/email.sh $inputfile "not exist" run_single_crest.sh $run_info
			$script_path/wait.sh $inputfile.fix.log 
			cat $inputfile | awk '$0 ~/^#/' > $basedir/struct/crest/$sample/vcf.header.$sam.crest
            cat $inputfile | awk '$0 !~ /^#/' >> $output/SV/$sam.crest.vcf
            cat $inputfile_filter | awk '$0 !~ /^#/' >> $output/SV/$sam.final.crest.vcf
            cat $input >> $basedir/struct/crest/$sample/$sam.crest
            cat $input_filter >> $basedir/struct/crest/$sample/$sam.final.crest
            rm $input
            rm $input_filter    
            rm $inputfile_filter
            rm $inputfile
        else
			cat $inputfile | awk '$0 ~/^#/' > $basedir/struct/crest/$sample/vcf.header.$sam.crest
            cat $inputfile | awk '$0 !~ /^#/' >> $output/SV/$sam.crest.vcf
            cat $inputfile_filter | awk '$0 !~ /^#/' >> $output/SV/$sam.final.crest.vcf
            cat $input >> $basedir/struct/crest/$sample/$sam.crest
            cat $input_filter >> $basedir/struct/crest/$sample/$sam.final.crest
            rm $input
            rm $input_filter    
            rm $inputfile_filter
            rm $inputfile
        fi
    done
    cat $basedir/struct/crest/$sample/vcf.header.$sam.crest $output/SV/$sam.crest.vcf > $output/SV/$sam.crest.vcf.temp
    mv $output/SV/$sam.crest.vcf.temp $output/SV/$sam.crest.vcf
    
    cat $basedir/struct/crest/$sample/vcf.header.$sam.crest $output/SV/$sam.final.crest.vcf > $output/SV/$sam.final.crest.vcf.temp
    mv $output/SV/$sam.final.crest.vcf.temp $output/SV/$sam.final.crest.vcf
	rm $basedir/struct/crest/$sample/vcf.header.$sam.crest

    if [ ! -s $output/SV/$sam.crest.vcf ]
    then
        $script_path/errorlog.sh $output/SV/$sam.crest.vcf summaryze_struct_single.sh ERROR "not created"
    	exit 1;
   	else
        file=`echo $inputargs | sed -e '/-V/s///g'`
        file=`echo $inputargs_filter | sed -e '/-V/s///g'`
    fi    
	echo `date`
fi    




