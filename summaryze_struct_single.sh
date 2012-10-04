#!/bin/sh

if [ $# != 3 ]
then	
	echo "Usage: <sample name> <base dir> <path to run_info file> ";
	exit 1
else
	set -x
	echo `date`
	sample=$1
	basedir=$2
	run_info=$3

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
		inputfile=$basedir/cnv/$sample/$sample.$chr.del.vcf
		input=$basedir/cnv/$sample/$sample.$chr.del.bed
		if [ ! -f $inputfile ]
		then
			touch $inputfile.fix.log
			$script_path/email.sh $inputfile "not exist" $JOB_NAME $JOB_ID $run_info
			$script_path/wait.sh $inputfile.fix.log 
			inputargs="-V $inputfile "$inputargs  
			cat $input >> $basedir/cnv/$sample/$sample.cnv.bed
			rm $input
		else
			inputargs="-V $inputfile "$inputargs  
			cat $input >> $basedir/cnv/$sample/$sample.cnv.bed
			rm $input
		fi
		
		inputfile=$basedir/cnv/$sample/$sample.$chr.dup.vcf 
		input=$basedir/cnv/$sample/$sample.$chr.dup.bed
		if [ ! -f $inputfile ]
		then	
			touch $inputfile.fix.log
			$script_path/email.sh $inputfile "not exist" $JOB_NAME $JOB_ID $run_info
			$script_path/wait.sh $inputfile.fix.log 
			inputargs="-V $inputfile "$inputargs  
			cat $input >> $basedir/cnv/$sample/$sample.cnv.bed
			rm $input
		else
			inputargs="-V $inputfile "$inputargs  
			cat $input >> $basedir/cnv/$sample/$sample.cnv.bed
			rm $input
		fi
		
		inputfile=$basedir/cnv/$sample/$sample.$chr.filter.del.vcf
		input=$basedir/cnv/$sample/$sample.$chr.filter.del.bed
		if [ ! -f $inputfile ]
		then	
			touch $inputfile.fix.log
			$script_path/email.sh $inputfile "not exist" $JOB_NAME $JOB_ID $run_info
			$script_path/wait.sh $inputfile.fix.log 
			cat $input >> $basedir/cnv/$sample/$sample.cnv.filter.bed
			rm $input
			inputargs_filter="-V $inputfile "$inputargs_filter  
		else
			cat $input >> $basedir/cnv/$sample/$sample.cnv.filter.bed
			rm $input
			inputargs_filter="-V $inputfile "$inputargs_filter  
		fi
		
		inputfile=$basedir/cnv/$sample/$sample.$chr.filter.dup.vcf
		input=$basedir/cnv/$sample/$sample.$chr.filter.dup.bed
		if [ ! -f $inputfile ]
		then
			touch $inputfile.fix.log
			$script_path/email.sh $inputfile "not exist" $JOB_NAME $JOB_ID $run_info
			$script_path/wait.sh $inputfile.fix.log 
			cat $input >> $basedir/cnv/$sample/$sample.cnv.filter.bed
			rm $input
			inputargs_filter="-V $inputfile "$inputargs_filter 
		else
			cat $input >> $basedir/cnv/$sample/$sample.cnv.filter.bed
			rm $input
			inputargs_filter="-V $inputfile "$inputargs_filter 
		fi
	done

	$script_path/combinevcf.sh "$inputargs" $output/SV/$sample.cnv.vcf $run_info yes
	$script_path/combinevcf.sh "$inputargs_filter" $output/SV/$sample.cnv.filter.vcf $run_info yes

	
    #Summaryzing Breakdancer
	inputargs=""
	input=""
	for chr in $chrs
	do
		inputfile=$basedir/struct/break/$sample/$sample.$chr.break.vcf 
		input=$basedir/struct/break/$sample/$sample.$chr.break
		if [ ! -s $inputfile ]
		then      
			touch $inputfile.fix.log
			$script_path/email.sh $inputfile "not exist" $JOB_NAME $JOB_ID $run_info
			$script_path/wait.sh $inputfile.fix.log 
			cat $inputfile | awk '$0 ~/^#/' > $basedir/struct/break/$sample/$sample.header.break
			cat $inputfile | awk '$0 !~ /^#/' >> $output/SV/$sample.break.vcf
			cat $input | awk '$0 !~ /^#/' >> $basedir/struct/break/$sample/$sample.break
			rm $input
			rm $inputfile
		else
			#inputargs="-V $inputfile "$inputargs
			cat $inputfile | awk '$0 ~/^#/' > $basedir/struct/break/$sample/$sample.header.break
			cat $inputfile | awk '$0 !~ /^#/' >> $output/SV/$sample.break.vcf
			cat $input | awk '$0 !~ /^#/' >> $basedir/struct/break/$sample/$sample.break
			rm $input
			rm $inputfile
		fi
	done
    cat $basedir/struct/break/$sample/$sample.inter.break | awk '$0 !~ /^#/' >> $basedir/struct/break/$sample/$sample.break
    cat $basedir/struct/break/$sample/$sample.inter.break.vcf |  awk '$0 !~ /^#/' >> $output/SV/$sample.break.vcf
	rm 	$basedir/struct/break/$sample/$sample.inter.break.vcf
	rm $basedir/struct/break/$sample/$sample.inter.break 
    cat $basedir/struct/break/$sample/$sample.header.break $output/SV/$sample.break.vcf > $output/SV/$sample.break.vcf.temp
	mv $output/SV/$sample.break.vcf.temp $output/SV/$sample.break.vcf
	rm $basedir/struct/break/$sample/$sample.header.break


	if [ ! -s $output/SV/$sample.break.vcf ]
	then
		$script_path/errorlog.sh $output/$sample.break.vcf summaryze_struct_single.sh ERROR "not created"
		exit 1;
	else
		file=`echo $inputargs | sed -e '/-V/s///g'`
	fi    
	#Summaryzing Crest

    inputargs=""
    inputargs_filter=""
    input=""
    input_filter=""
    for chr in $chrs
    do
        inputfile=$basedir/struct/crest/$sample/$sample.$chr.raw.vcf
        inputfile_filter=$basedir/struct/crest/$sample/$sample.$chr.filter.vcf
        input=$basedir/struct/crest/$sample/$sample.$chr.predSV.txt
        input_filter=$basedir/struct/crest/$sample/$sample.$chr.filter.predSV.txt
        if [ ! -s $inputfile ]
        then      
            touch $inputfile.fix.log
			$script_path/email.sh $inputfile "not exist" $JOB_NAME $JOB_ID $run_info
			$script_path/wait.sh $inputfile.fix.log 
			cat $inputfile | awk '$0 ~/^#/' > $basedir/struct/crest/$sample/vcf.header.$sample.crest
            cat $inputfile | awk '$0 !~ /^#/' >> $output/SV/$sample.raw.crest.vcf
            cat $inputfile_filter | awk '$0 !~ /^#/' >> $output/SV/$sample.filter.crest.vcf
            cat $input >> $basedir/struct/crest/$sample/$sample.raw.crest
            cat $input_filter >> $basedir/struct/crest/$sample/$sample.filter.crest
            rm $input
            rm $input_filter    
            rm $inputfile_filter
            rm $inputfile
        else
			cat $inputfile | awk '$0 ~/^#/' > $basedir/struct/crest/$sample/vcf.header.$sample.crest
            cat $inputfile | awk '$0 !~ /^#/' >> $output/SV/$sample.raw.crest.vcf
            cat $inputfile_filter | awk '$0 !~ /^#/' >> $output/SV/$sample.filter.crest.vcf
            cat $input >> $basedir/struct/crest/$sample/$sample.raw.crest
            cat $input_filter >> $basedir/struct/crest/$sample/$sample.filter.crest
            rm $input
            rm $input_filter    
            rm $inputfile_filter
            rm $inputfile
        fi
    done
    cat $basedir/struct/crest/$sample/vcf.header.$sample.crest $output/SV/$sample.raw.crest.vcf > $output/SV/$sample.raw.crest.vcf.temp
    mv $output/SV/$sample.raw.crest.vcf.temp $output/SV/$sample.raw.crest.vcf
    
    cat $basedir/struct/crest/$sample/vcf.header.$sample.crest $output/SV/$sample.filter.crest.vcf > $output/SV/$sample.filter.crest.vcf.temp
    mv $output/SV/$sample.filter.crest.vcf.temp $output/SV/$sample.filter.crest.vcf
	rm $basedir/struct/crest/$sample/vcf.header.$sample.crest

    if [ ! -s $output/SV/$sample.raw.crest.vcf ]
    then
        $script_path/errorlog.sh $output/SV/$sample.raw.crest.vcf summaryze_struct_single.sh ERROR "not created"
    	exit 1;
   	else
        file=`echo $inputargs | sed -e '/-V/s///g'`
        file=`echo $inputargs_filter | sed -e '/-V/s///g'`
    fi    
	echo `date`
fi    




