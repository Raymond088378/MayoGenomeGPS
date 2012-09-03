#!/bin/bash

if [ $# != 4 ]
then
    echo -e "Usage: combine vcfs <input files><output vcf ><run info><to delete input files(yes/no)"
else
    set -x
    echo `date`
    input=$1
    output=$2
    run_info=$3
    flag=`echo $4 | tr "[a-z]" "[A-Z]"`
    
    tool_info=$( cat $run_info | grep -w '^TOOL_INFO' | cut -d '=' -f2)
	script_path=$( cat $tool_info | grep -w '^WHOLEGENOME_PATH' | cut -d '=' -f2)
    vcftools=$( cat $tool_info | grep -w '^VCFTOOLS' | cut -d '=' -f2)
	perllib=$( cat $tool_info | grep -w '^PERLLIB_VCF' | cut -d '=' -f2)
	tabix=$( cat $tool_info | grep -w '^TABIX' | cut -d '=' -f2)
	ref=$( cat $tool_info | grep -w '^REF_GENOME' | cut -d '=' -f2)
	export PERL5LIB=$PERL5LIB:$perllib
	export PATH=$tabix/:$PATH
	
	args=""
	for file in $input
	do
		if [ ! -s $file.gz ]
		then
			$tabix/bgzip $file
		fi
		$tabix/tabix -p vcf $file.gz
		args=$args" $file.gz"
	done	
	
	$vcftools/bin/vcf-concat $args > $output
    
    if [ ! -s $output ]
    then
        $script_path/errorlog.sh $output concatvcf.sh ERROR "failed to create"
    else
        perl $script_path/vcfsort.pl $ref.fai $output > $output.temp
        mv $output.temp $output
		for file in $input
		do
			$tabix/bgzip -d $file.gz
			rm $file.gz.tbi
		done	
        if [ $flag == "YES" ]
        then
            rm $input
            sample=`echo $input | tr " " "\n" | awk '{print $0".idx"}'`
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
