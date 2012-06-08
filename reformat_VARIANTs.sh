#!/bin/sh

#	INFO
#	reformat the inputs to chop into chromsome to make it faster fro vraiant module

if [ $# != 4 ]
then
    echo "usage: <output><sample><run_info><marker>";
else	
    set -x
    echo `date`
    output=$1
    sample=$2
    run_info=$3
    marker=$4
    
    input=$( cat $run_info | grep -w '^INPUT_DIR' | cut -d '=' -f2)
    tool_info=$( cat $run_info | grep -w '^TOOL_INFO' | cut -d '=' -f2)
    sample_info=$( cat $run_info | grep -w '^SAMPLE_INFO' | cut -d '=' -f2)
    variant_type=$( cat $run_info | grep -w '^VARIANT_TYPE' | cut -d '=' -f2)
    script_path=$( cat $tool_info | grep -w '^WHOLEGENOME_PATH' | cut -d '=' -f2)
    dbsnp_rsids_snv=$( cat $tool_info | grep -w '^dbSNP_SNV_rsIDs' | cut -d '=' -f2)
    SNV_caller=$( cat $run_info | grep -w '^SNV_CALLER' | cut -d '=' -f2)
    variant_type=`echo "$variant_type" | tr "[a-z]" "[A-Z]"`
    dbsnp_rsids_indel=$( cat $tool_info | grep -w '^dbSNP_INDEL_rsIDs' | cut -d '=' -f2)
    chrs=$( cat $run_info | grep -w '^CHRINDEX' | cut -d '=' -f2 | tr ":" " " )
    sample_info=$( cat $run_info | grep -w '^SAMPLE_INFO' | cut -d '=' -f2)
    distance=$( cat $tool_info | grep -w '^SNP_DISTANCE_INDEL' | cut -d '=' -f2 )
	blat=$( cat $tool_info | grep -w '^BLAT' | cut -d '=' -f2 )
    blat_port=$( cat $tool_info | grep -w '^BLAT_PORT' | cut -d '=' -f2 )
    blat_ref=$( cat $tool_info | grep -w '^BLAT_REF' | cut -d '=' -f2 )
    blat_server=$( cat $tool_info | grep -w '^BLAT_SERVER' | cut -d '=' -f2 )
    window_blat=$( cat $tool_info | grep -w '^WINDOW_BLAT' | cut -d '=' -f2 )
	javahome=$( cat $tool_info | grep -w '^JAVA_HOME' | cut -d '=' -f2 )
	threads=$( cat $tool_info | grep -w '^THREADS' | cut -d '=' -f2 )
	samtools=$( cat $tool_info | grep -w '^SAMTOOLS' | cut -d '=' -f2 )
	ref=$( cat $tool_info | grep -w '^REF_GENOME' | cut -d '=' -f2)
	perllib=$( cat $tool_info | grep -w '^PERLLIB' | cut -d '=' -f2)
	export PERL5LIB=$perllib:$PERL5LIB
	export PATH=$PERL5LIB:$PATH
	range=20000
    let blat_port+=$RANDOM%range
    status=`$blat/gfServer status $blat_server $blat_port | wc -l`;
    if [ "$status" -le 1 ]
    then
		$blat/gfServer start $blat_server $blat_port -log=$output/$sample.blat.log $blat_ref  &
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
            rm $output/$sample.blat.log
            $blat/gfServer start $blat_server $blat_port -log=$output/$sample.blat.log $blat_ref  &
            sleep 3m
        fi
		status=`$blat/gfServer status $blat_server $blat_port | wc -l`;
    done
	
	
    if [ $marker -eq 2 ]
    then
        snv_file=$( cat $sample_info | grep -w SNV:${sample} | cut -d '=' -f2)
        indel_file=$( cat $sample_info | grep -w INDEL:${sample} | cut -d '=' -f2)
        
		if [ "$input/$snv_file" != "$input/$indel_file" ]
		then
			## format the vcf file to text delimited file
			## determine if the file is vcf or text input
			type=`cat $input/$snv_file | head -1 | awk '{if ($0 ~ /^##/) print "vcf";else print "txt"}'`
			if [ $type == "txt" ]
			then
				perl $script_path/convert_txt_vcf.pl $input/$snv_file $sample > $output/$sample.SNV.vcf
				perl $script_path/convert_txt_vcf.pl $input/$indel_file $sample > $output/$sample.INDEL.vcf
			else
				cp $input/$snv_file $output/$sample.SNV.vcf
				cp $input/$indel_file $output/$sample.INDEL.vcf 
			fi
			n=`cat $output/$sample.SNV.vcf |  awk '$0 ~ /^##INFO=<ID=ED/' | wc -l`
			if [ $n == 0 ]
			then
				$script_path/vcf_blat_verify.pl -i $output/$sample.SNV.vcf -o $output/$sample.SNV.vcf.tmp -w $window_blat -b $blat -r $ref -sam $samtools -br $blat_ref -bs $blat_server -bp $blat_port -th $threads
				perl $script_path/vcfsort.pl $ref.fai $output/$sample.SNV.vcf.tmp > $output/$sample.SNV.vcf
				rm $output/$sample.SNV.vcf.tmp
			else
				perl $script_path/vcfsort.pl $ref.fai $output/$sample.SNV.vcf > $output/$sample.SNV.vcf.tmp
				mv $output/$sample.SNV.vcf.tmp $output/$sample.SNV.vcf
			fi	
			n=`cat $output/$sample.INDEL.vcf |  awk '$0 ~ /^##INFO=<ID=ED/' | wc -l`
			if [ $n == 0 ]
			then
				$script_path/vcf_blat_verify.pl -i $output/$sample.INDEL.vcf -o $output/$sample.INDEL.vcf.tmp -w $window_blat -b $blat -r $ref -sam $samtools -br $blat_ref -bs $blat_server -bp $blat_port -th $threads
				perl $script_path/vcfsort.pl $ref.fai $output/$sample.INDEL.vcf.tmp > $output/$sample.INDEL.vcf
				rm $output/$sample.INDEL.vcf.tmp 
			else
				perl $script_path/vcfsort.pl $ref.fai $output/$sample.INDEL.vcf > $output/$sample.INDEL.vcf.tmp
				mv $output/$sample.INDEL.vcf.tmp $output/$sample.INDEL.vcf 
			fi
		else
			type=`cat $input/$snv_file | head -1 | awk '{if ($0 ~ /^##/) print "vcf";else print "txt"}'`
			if [ $type == "txt" ]
			then
				perl $script_path/convert_txt_vcf.pl $input/$snv_file $sample > $output/$sample.vcf
			else
				cp $input/$snv_file $output/$sample.vcf
			fi
			n=`cat $output/$sample.vcf |  awk '$0 ~ /^##INFO=<ID=ED/' | wc -l`
			if [ $n == 0 ]
			then
				$script_path/vcf_blat_verify.pl -i $output/$sample.vcf -o $output/$sample.vcf.tmp -w $window_blat -b $blat -r $ref -sam $samtools -br $blat_ref -bs $blat_server -bp $blat_port -th $threads
				perl $script_path/vcfsort.pl $ref.fai $output/$sample.vcf.tmp > $output/$sample.vcf
				rm $output/$sample.vcf.tmp
			else
				perl $script_path/vcfsort.pl $ref.fai $output/$sample.vcf > $output/$sample.vcf.tmp
				mv $output/$sample.vcf.tmp $output/$sample.vcf
			fi
			perl $script_path/vcf_to_variant_vcf.pl -i $output/$sample.vcf -v $output/$sample.SNV.vcf -l $output/$sample.INDEL.vcf
			rm $output/$sample.vcf
		fi	
			
		for chr in $chrs
        do
            perl $script_path/vcf_to_variant_vcf.pl -i $output/$sample.SNV.vcf -v $output/$sample.variants.chr$chr.SNV.filter.i.c.vcf -l $output/$sample.variants.chr$chr.INDEL.filter.i.c.vcf -c chr$chr -s $sample
            rm $output/$sample.variants.chr$chr.INDEL.filter.i.c.vcf
            perl $script_path/vcf_to_variant_vcf.pl -i $output/$sample.INDEL.vcf -v $output/$sample.variants.chr$chr.SNV.filter.i.c.vcf.temp -l $output/$sample.variants.chr$chr.INDEL.filter.i.c.vcf -c chr$chr -s $sample
            rm $output/$sample.variants.chr$chr.SNV.filter.i.c.vcf.temp
            cat $output/$sample.variants.chr$chr.INDEL.filter.i.c.vcf | awk 'BEGIN {OFS="\t"} {if ($0 ~ /^#/) print $0; else print $1,$2,$3,$4,$5,$6,$7,$8";CAPTURE=1",$9,$10;}' > $output/$sample.variants.chr$chr.INDEL.filter.i.c.vcf.tmp
            perl $script_path/add_format_field_vcf.pl $output/$sample.variants.chr$chr.INDEL.filter.i.c.vcf.tmp INDEL > $output/$sample.variants.chr$chr.INDEL.filter.i.c.vcf
            rm $output/$sample.variants.chr$chr.INDEL.filter.i.c.vcf.tmp
            if [ `cat $output/$sample.variants.chr$chr.INDEL.filter.i.c.vcf | awk '$0 !~ /^#/' | wc -l` -ge 1 ]
			then
				perl $script_path/markSnv_IndelnPos.pl -s $output/$sample.variants.chr$chr.SNV.filter.i.c.vcf -i $output/$sample.variants.chr$chr.INDEL.filter.i.c.vcf -n $distance -o $output/$sample.variants.chr$chr.SNV.filter.i.c.pos.vcf
				cat $output/$sample.variants.chr$chr.SNV.filter.i.c.pos.vcf | awk 'BEGIN {OFS="\t"} {if ($0 ~ /^#/) print $0; else print $1,$2,$3,$4,$5,$6,$7,$8";CAPTURE=1;CLOSE2INDEL="$NF,$9,$10;}' > $output/$sample.variants.chr$chr.SNV.filter.i.c.vcf  
			else
				cat $output/$sample.variants.chr$chr.SNV.filter.i.c.vcf | awk 'BEGIN {OFS="\t"} {if ($0 ~ /^#/) print $0; else print $1,$2,$3,$4,$5,$6,$7,$8";CAPTURE=1;CLOSE2INDEL=0",$9,$10;}' > $output/$sample.variants.chr$chr.SNV.filter.i.c.vcf.tmp
				mv $output/$sample.variants.chr$chr.SNV.filter.i.c.vcf.tmp $output/$sample.variants.chr$chr.SNV.filter.i.c.vcf	
			fi	
			perl $script_path/add_format_field_vcf.pl $output/$sample.variants.chr$chr.SNV.filter.i.c.vcf SNV > $output/$sample.variants.chr$chr.SNV.filter.i.c.pos.vcf 
            mv $output/$sample.variants.chr$chr.SNV.filter.i.c.pos.vcf $output/$sample.variants.chr$chr.SNV.filter.i.c.vcf		
        done
        rm $output/$sample.SNV.vcf $output/$sample.INDEL.vcf 
    else
        if [ $variant_type == "SNV" ]
        then
            snv_file=$( cat $sample_info | grep -w SNV:${sample} | cut -d '=' -f2)
            type=`cat $input/$snv_file | head -1 | awk '{if ($0 ~ /^##/) print "vcf";else print "txt"}'`
            if [ $type == "txt" ]
            then
                perl $script_path/convert_txt_vcf.pl $input/$snv_file $sample > $output/$sample.SNV.vcf
            else
                cp $input/$snv_file $output/$sample.SNV.vcf
            fi 
			n=`cat $output/$sample.SNV.vcf |  awk '$0 ~ /^##INFO=<ID=ED/' | wc -l`
			if [ $n == 0 ]
			then
				$script_path/vcf_blat_verify.pl -i $output/$sample.SNV.vcf -o $output/$sample.SNV.vcf.tmp -w $window_blat -b $blat -r $ref -sam $samtools -br $blat_ref -bs $blat_server -bp $blat_port -th $threads
				perl $script_path/vcfsort.pl $ref.fai $output/$sample.SNV.vcf.tmp > $output/$sample.SNV.vcf
				rm $output/$sample.SNV.vcf.tmp
			else
				perl $script_path/vcfsort.pl $ref.fai $output/$sample.SNV.vcf > $output/$sample.SNV.vcf.tmp
				mv $output/$sample.SNV.vcf.tmp $output/$sample.SNV.vcf
			fi	
            for chr in $chrs	
            do
                perl $script_path/vcf_to_variant_vcf.pl -i $output/$sample.SNV.vcf -v $output/$sample.variants.chr$chr.SNV.filter.i.c.vcf -l $output/$sample.variants.chr$chr.INDEL.filter.i.c.vcf -c chr$chr -s $sample
                rm $output/$sample.variants.chr$chr.INDEL.filter.i.c.vcf
                cat $output/$sample.variants.chr$chr.SNV.filter.i.c.vcf | awk 'BEGIN {OFS="\t"} {if ($0 ~ /^#/) print $0; else print $1,$2,$3,$4,$5,$6,$7,$8";CAPTURE=1;CLOSE2INDEL=0",$9,$10;}' > $output/$sample.variants.chr$chr.SNV.filter.i.c.pos.vcf
				perl $script_path/add_format_field_vcf.pl $output/$sample.variants.chr$chr.SNV.filter.i.c.pos.vcf SNV > $output/$sample.variants.chr$chr.SNV.filter.i.c.vcf
                rm $output/$sample.variants.chr$chr.SNV.filter.i.c.pos.vcf			
            done
            rm $output/$sample.SNV.vcf
        elif [ $variant_type == "INDEL" ]
        then
            indel_file=$( cat $sample_info | grep -w INDEL:${sample} | cut -d '=' -f2)
            type=`cat $input/$indel_file | head -1 | awk '{if ($0 ~ /^##/) print "vcf";else print "txt"}'`
            if [ $type == "txt" ]
            then
                perl $script_path/convert_txt_vcf.pl $input/$indel_file $sample > $output/$sample.INDEL.vcf 
            else
                cp $input/$indel_file $output/$sample.INDEL.vcf 
            fi    
            n=`cat $output/$sample.INDEL.vcf |  awk '$0 ~ /^##INFO=<ID=ED/' | wc -l`
			if [ $n == 0 ]
			then
				$script_path/vcf_blat_verify.pl -i $output/$sample.INDEL.vcf -o $output/$sample.INDEL.vcf.tmp -w $window_blat -b $blat -r $ref -sam $samtools -br $blat_ref -bs $blat_server -bp $blat_port -th $threads
				perl $script_path/vcfsort.pl $ref.fai $output/$sample.INDEL.vcf.tmp > $output/$sample.INDEL.vcf
				rm $output/$sample.INDEL.vcf.tmp 
			else
				perl $script_path/vcfsort.pl $ref.fai $output/$sample.INDEL.vcf > $output/$sample.INDEL.vcf.tmp
				mv $output/$sample.INDEL.vcf.tmp $output/$sample.INDEL.vcf 
			fi
			for chr in $chrs		
            do
                perl $script_path/vcf_to_variant_vcf.pl -i $output/$sample.INDEL.vcf -v $output/$sample.variants.chr$chr.SNV.filter.i.c.vcf -l $output/$sample.variants.chr$chr.INDEL.filter.i.c.vcf -c chr$chr -s $sample
                rm $output/$sample.variants.chr$chr.SNV.filter.i.c.vcf
                cat $output/$sample.variants.chr$chr.INDEL.filter.i.c.vcf | awk 'BEGIN {OFS="\t"} {if ($0 ~ /^#/) print $0; else print $1,$2,$3,$4,$5,$6,$7,$8";CAPTURE=1",$9,$10;}' > $output/$sample.variants.chr$chr.INDEL.filter.i.c.vcf.tmp
				perl $script_path/add_format_field_vcf.pl $output/$sample.variants.chr$chr.INDEL.filter.i.c.vcf.tmp INDEL > $output/$sample.variants.chr$chr.INDEL.filter.i.c.vcf  
                rm $output/$sample.variants.chr$chr.INDEL.filter.i.c.vcf.tmp
            done
            rm $output/$sample.INDEL.vcf 
        fi		
    fi
    rm $output/$sample.blat.log
    echo `date`	
fi	
	
	
	
	
	
	
	
		
