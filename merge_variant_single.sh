#!/bin/sh

########################################################
###### 	Merges variants from vcf files by chromosome

######		Program:			merge_variant_group.sh
######		Date:				12/13/2011
######		Summary:			Using PICARD to sort and mark duplicates in bam 
######		Input files:		$1	=	/path/to/input directory
######					$2	=	sample name
######					$3	=	/path/to/run_info.txt
########################################################

if [ $# != 4 ];
then
    echo "Usage: </path/to/input directory> <sample name> </path/to/run_info.txt>";
else
    set -x
    echo `date`
    input=$1
    sample=$2
    out=$3
    run_info=$4
    
########################################################	
######		Reading run_info.txt and assigning to variables
    tool_info=$( cat $run_info | grep -w '^TOOL_INFO' | cut -d '=' -f2)
    picard=$( cat $tool_info | grep -w '^PICARD' | cut -d '=' -f2 ) 
    java=$( cat $tool_info | grep -w '^JAVA' | cut -d '=' -f2)
    chrs=$( cat $run_info | grep -w '^CHRINDEX' | cut -d '=' -f2 | tr ":" "\n" )
    ref=$( cat $tool_info | grep -w '^REF_GENOME' | cut -d '=' -f2)
    gatk=$( cat $tool_info | grep -w '^GATK' | cut -d '=' -f2)
    script_path=$( cat $tool_info | grep -w '^WHOLEGENOME_PATH' | cut -d '=' -f2 )

    output=$( cat $run_info | grep -w '^BASE_OUTPUT_DIR' | cut -d '=' -f2)
    PI=$( cat $run_info | grep -w '^PI' | cut -d '=' -f2)
    tool=$( cat $run_info | grep -w '^TYPE' | cut -d '=' -f2 | tr "[A-Z]" "[a-z]" )
    run_num=$( cat $run_info | grep -w '^OUTPUT_FOLDER' | cut -d '=' -f2)
    filter_variants=$( cat $tool_info | grep -w '^VARIANT_FILTER' | cut -d '=' -f2 | tr "[a-z]" "[A-Z]")
    blat=$( cat $tool_info | grep -w '^BLAT' | cut -d '=' -f2 )
    blat_port=$( cat $tool_info | grep -w '^BLAT_PORT' | cut -d '=' -f2 )
    blat_ref=$( cat $tool_info | grep -w '^BLAT_REF' | cut -d '=' -f2 )
    blat_server=$( cat $tool_info | grep -w '^BLAT_SERVER' | cut -d '=' -f2 )
    window_blat=$( cat $tool_info | grep -w '^WINDOW_BLAT' | cut -d '=' -f2 )

    range=20000
    let blat_port+=$RANDOM%range
    status=`$blat/gfServer status localhost $blat_port | wc -l`;
    if [ "$status" -le 1 ]
    then
	$blat/gfServer start localhost $blat_port -log=$out/$sample.variants.raw.vcf.blat.log $blat_ref  &
	sleep 2m
    fi
    status=`$blat/gfServer status localhost $blat_port | wc -l`;

    if [ "$status" -eq 0 ]
    then
        blat_port=$( cat $tool_info | grep -w '^BLAT_PORT' | cut -d '=' -f2 )
        range=20000
        let blat_port+=$RANDOM%range
        status=`$blat/gfServer status localhost $blat_port | wc -l`;
        if [ "$status" -le 1 ]
        then
            rm $out/$sample.variants.raw.vcf.blat.log
            $blat/gfServer start localhost $blat_port -log=$out/$sample.variants.raw.vcf.blat.log $blat_ref  &
            sleep 2m
        fi
    fi
    
    inputargs=""
    for i in $chrs
    do
      inputfile=$input/$sample/$sample.variants.chr$i.raw.vcf 
      if [ ! -s $inputfile ]
      then	
          echo "ERROR :merge_variant_single.sh Variant file for sample $sample, chromosome $i: $inputfile does not exist "
          exit 1;
      else
          inputargs="-V $inputfile "$inputargs
      fi
    done

    $script_path/combinevcf.sh "$inputargs" $out/$sample.variants.raw.vcf $run_info no
    perl $script_path/vcf_blat_verify.pl -i $out/$sample.variants.raw.vcf -o $out/$sample.variants.raw.vcf.tmp -w $window_blat -b $blat -r $ref -br $blat_ref -bs $blat_server -bp $blat_port
    mv $out/$sample.variants.raw.vcf.tmp $out/$sample.variants.raw.vcf
    
    ### filter the variant calls
    if [ $filter_variants == "YES" ]
    then
        $script_path/filter_variant_vqsr.sh $out/$sample.variants.raw.vcf $out/$sample.variants.filter.vcf BOTH $run_info
    else
        cp $out/$sample.variants.raw.vcf $out/$sample.variants.filter.vcf
    fi
    
    if [ ! -s $out/$sample.variants.filter.vcf ]
    then
        echo "ERROR: failed to generate the filterd vcf for $sample"
        exit 1;
    else
        for i in $chrs
        do
            cat $out/$sample.variants.filter.vcf | awk -v num=chr${i} '$0 ~ /^#/ || $1 == num' > $input/$sample/$sample.variants.chr$i.filter.vcf 
        done
    fi  
	echo `date`	
fi  
