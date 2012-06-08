#!/bin/sh

if [ $# != 5 ]
then
    echo "usage: <output><sample><run_info><marker>";
else	
    set -x
    echo `date`
    output=$1
    reports=$2
    sample=$3
    run_info=$4
    marker=$5
    
    input=$( cat $run_info | grep -w '^INPUT_DIR' | cut -d '=' -f2)
    tool_info=$( cat $run_info | grep -w '^TOOL_INFO' | cut -d '=' -f2)
    sample_info=$( cat $run_info | grep -w '^SAMPLE_INFO' | cut -d '=' -f2)
    variant_type=$( cat $run_info | grep -w '^VARIANT_TYPE' | cut -d '=' -f2| tr "[a-z]" "[A-Z]")
    script_path=$( cat $tool_info | grep -w '^WHOLEGENOME_PATH' | cut -d '=' -f2)
    dbsnp_rsids_indel=$( cat $tool_info | grep -w '^dbSNP_INDEL_rsIDs' | cut -d '=' -f2)
    chrs=$( cat $run_info | grep -w '^CHRINDEX' | cut -d '=' -f2 | tr ":" " " )
    threads=$( cat $tool_info | grep -w '^THREADS' | cut -d '=' -f2 )
    samtools=$( cat $tool_info | grep -w '^SAMTOOLS' | cut -d '=' -f2 )
    vcftools=$( cat $tool_info | grep -w '^VCFTOOLS' | cut -d '=' -f2)
    perllib=$( cat $tool_info | grep -w '^PERLLIB_VCF' | cut -d '=' -f2)
    blat=$( cat $tool_info | grep -w '^BLAT' | cut -d '=' -f2 )
    blat_port=$( cat $tool_info | grep -w '^BLAT_PORT' | cut -d '=' -f2 )
    blat_ref=$( cat $tool_info | grep -w '^BLAT_REF' | cut -d '=' -f2 )
    blat_server=$( cat $tool_info | grep -w '^BLAT_SERVER' | cut -d '=' -f2 )
    window_blat=$( cat $tool_info | grep -w '^WINDOW_BLAT' | cut -d '=' -f2 )
    threads=$( cat $tool_info | grep -w '^THREADS' | cut -d '=' -f2 )
    ref=$( cat $tool_info | grep -w '^REF_GENOME' | cut -d '=' -f2)
    perllib=$( cat $tool_info | grep -w '^PERLLIB' | cut -d '=' -f2)
	export PERL5LIB=$perllib:$PERL5LIB
	export PATH=$PERL5LIB:$PATH
	
	mkdir -p $output/$sample/
    
    range=20000
    let blat_port+=$RANDOM%range
    status=`$blat/gfServer status localhost $blat_port | wc -l`;
    if [ "$status" -le 1 ]
    then
		$blat/gfServer start $blat_server $blat_port -log=$output/$sample/$sample.blat.log $blat_ref  &
		sleep 3m
    fi
    status=`$blat/gfServer status $blat_server $blat_port | wc -l`;

    while [ "$status" -le 1 ]
    do
        blat_port=$( cat $tool_info | grep -w '^BLAT_PORT' | cut -d '=' -f2 )
        range=20000
        let blat_port+=$RANDOM%range
        status=`$blat/gfServer status $blat_server $blat_port | wc -l`;
        if [ "$status" -le 1 ]
        then
            if [ -s $output/$sample/$sample.blat.log ]
			then
				rm $output/$sample/$sample.blat.log
            fi
			$blat/gfServer start $blat_server $blat_port -log=$output/$sample/$sample.blat.log $blat_ref  &
            sleep 3m
        fi
		status=`$blat/gfServer status $blat_server $blat_port | wc -l`;
    done
    
    
    if [ $marker == 2 ]
    then
        snv_file=$( cat $sample_info | grep -w SNV:${sample} | cut -d '=' -f2)
        indel_file=$( cat $sample_info | grep -w INDEL:${sample} | cut -d '=' -f2)	
	
        if [ "$input/$snv_file" == "$input/$indel_file" ]
        then
			n=`cat $input/$snv_file |  awk '$0 ~ /^##INFO=<ID=ED/' | wc -l`
			if [ $n == 0 ]
			then
				$script_path/vcf_blat_verify.pl -i $input/$snv_file -o $reports/$sample.variants.raw.vcf -w $window_blat -b $blat -r $ref -br $blat_ref -bs $blat_server -bp $blat_port -sam $samtools -th $threads
			else
				cp $input/$snv_file $reports/$sample.variants.raw.vcf
			fi
			cp $reports/$sample.variants.raw.vcf $reports/$sample.variants.filter.vcf
			for chr in $chrs
			do
				cat $reports/$sample.variants.filter.vcf | awk -v num=chr${chr} '$0 ~ /^#/ || $1 == num' > $output/$sample/$sample.variants.chr$chr.filter.vcf
			done
		else
			### concat both the vcf files to get one vcf file
			snv=`cat $input/$snv_file | awk '$0 ~ /^##INFO=<ID=ED/' | wc -l`
			if [ $snv == 0 ]
			then
				$script_path/vcf_blat_verify.pl -i $input/$snv_file -o $reports/$sample.variants.SNV.raw.vcf -w $window_blat -b $blat -r $ref -br $blat_ref -bs $blat_server -bp $blat_port -sam $samtools -th $threads
			else
				cp $input/$snv_file $reports/$sample.variants.SNV.raw.vcf
			fi
			indel=`cat $input/$indel_file | awk '$0 ~ /^##INFO=<ID=ED/' | wc -l`
			if [ $indel == 0 ]
			then
				$script_path/vcf_blat_verify.pl -i $input/$indel_file -o $reports/$sample.variants.INDEL.raw.vcf -w $window_blat -b $blat -r $ref -br $blat_ref -bs $blat_server -bp $blat_port -sam $samtools -th $threads
			else
				cp $input/$indel_file $reports/$sample.variants.INDEL.raw.vcf
			fi	
				in="-V $reports/$sample.variants.SNV.raw.vcf $reports/$sample.variants.INDEL.raw.vcf"
				$script_path/combinevcf.sh "$in" $reports/$sample.variants.raw.vcf $run_info yes
				cp $reports/$sample.variants.raw.vcf $reports/$sample.variants.filter.vcf
			for chr in $chrs
            do
               cat $reports/$sample.variants.filter.vcf | awk -v num=chr${chr} '$0 ~ /^#/ || $1 == num' > $output/$sample/$sample.variants.chr$chr.filter.vcf 
            done    
        fi
    else
        if [ $variant_type == "SNV" ]
        then
            snv_file=$( cat $sample_info | grep -w SNV:${sample} | cut -d '=' -f2)
            snv=`cat $input/$snv_file | awk '$0 ~ /^##INFO=<ID=ED/' | wc -l`
	    if [ $snv == 0 ]
	    then
		$script_path/vcf_blat_verify.pl -i $input/$snv_file -o $reports/$sample.variants.raw.vcf -w $window_blat -b $blat -r $ref -br $blat_ref -bs $blat_server -bp $blat_port -sam $samtools -th $threads
	    else
		cp $input/$snv_file $reports/$sample.variants.raw.vcf
	    fi
            cp $reports/$sample.variants.raw.vcf $reports/$sample.variants.filter.vcf
            for chr in $chrs
            do
                cat $reports/$sample.variants.filter.vcf | awk -v num=chr${chr} '$0 ~ /^#/ || $1 == num' > $output/$sample/$sample.variants.chr$chr.filter.vcf 
            done 
        elif [ $variant_type == "INDEL" ]
        then
            indel_file=$( cat $sample_info | grep -w INDEL:${sample} | cut -d '=' -f2)
            indel=`cat $input/$indel_file | awk '$0 ~ /^##INFO=<ID=ED/' | wc -l`
	    if [ $indel == 0 ]
	    then
		$script_path/vcf_blat_verify.pl -i $input/$indel_file -o $reports/$sample.variants.raw.vcf -w $window_blat -b $blat -r $ref -br $blat_ref -bs $blat_server -bp $blat_port -sam $samtools -th $threads
	    else
		cp $input/$indel_file $reports/$sample.variants.raw.vcf
	    fi
            cp $reports/$sample.variants.raw.vcf $reports/$sample.variants.filter.vcf
            for chr in $chrs
            do
                cat $reports/$sample.variants.filter.vcf | awk -v num=chr${chr} '$0 ~ /^#/ || $1 == num' > $output/$sample/$sample.variants.chr$chr.filter.vcf 
            done 
        fi
    fi
    echo `date`
fi    
            