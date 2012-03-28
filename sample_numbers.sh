#!/bin/sh
## script to get statistical numbers for each sample or pair

if [ $# != 2 ]
then
    echo -e "Usage: scrip to get statistical numbers for a sample\n <input dir> <run info>"
else 
    set -x
    echo `date`
    input_dir=$1
    run_info=$2
#SGE_TASK_ID=1
##############################################################		
    tool_info=$( cat $run_info | grep -w '^TOOL_INFO' | cut -d '=' -f2)
	sample_info=$( cat $run_info | grep -w '^SAMPLE_INFO' | cut -d '=' -f2)
    sample=$(cat $run_info | grep -w '^SAMPLENAMES' | cut -d '=' -f2 | tr ":" "\n" | head -n $SGE_TASK_ID | tail -n 1)
    group=$(cat $run_info | grep -w '^GROUPNAMES' | cut -d '=' -f2 | tr ":" "\n" | head -n $SGE_TASK_ID | tail -n 1)
    aligner=$( cat $run_info | grep -w '^ALIGNER' | cut -d '=' -f2)
	variant_type=$( cat $run_info | grep -w '^VARIANT_TYPE' | cut -d '=' -f2)
    run_num=$( cat $run_info | grep -w '^OUTPUT_FOLDER' | cut -d '=' -f2)
    flowcell=`echo $run_num | awk -F'_' '{print $NF}' | sed 's/.\(.*\)/\1/'`
	markdup=$( cat $run_info | grep -w '^MARKDUP' | cut -d '=' -f2| tr "[A-Z]" "[a-z]")
    chrs=$( cat $run_info | grep -w '^CHRINDEX' | cut -d '=' -f2 | tr ":" " ")
    caller=$( cat $run_info | grep -w '^SNV_CALLER' | cut -d '=' -f2)
    multi_sample=$( cat $run_info | grep -w '^MULTISAMPLE' | cut -d '=' -f2)  
    script_path=$( cat $tool_info | grep -w '^WHOLEGENOME_PATH' | cut -d '=' -f2 )
	tool=$( cat $run_info | grep -w '^TYPE' | cut -d '=' -f2| tr "[A-Z]" "[a-z]")
	analysis=$( cat $run_info | grep -w '^ANALYSIS' | cut -d '=' -f2| tr "[A-Z]" "[a-z]")
	version=$( cat $run_info | grep -w '^VERSION' | cut -d '=' -f2)
##############################################################		
##############################################################		
    
    numbers=$input_dir/numbers
	if [ $analysis == "mayo" -o $analysis == "realign-mayo" ]
	then
		pos=$( cat $run_info | grep -w '^SAMPLENAMES' | cut -d '=' -f2 | tr ":" "\n" | grep -n $sample | cut -d ":" -f1)
		lanes=$( cat $run_info | grep -w '^LANEINDEX' | cut -d '=' -f2 | tr ":" "\n" | head -n $pos | tail -n 1 | tr "," " ")
		i=1
		for lane in $lanes
		do
			index=$( cat $run_info | grep -w '^LABINDEXES' | cut -d '=' -f2 | tr ":" "\n" | head -n $pos | tail -n 1 | tr "," "\n" | head -n $i | tail -n 1)
			if [ $index == "-" ]
			then
				$java/java -jar $script_path/AddSecondaryAnalysis.jar -p $script_path/AddSecondaryAnalysis.properties -l $lane -f $flowcell -r $run_num -s Statistics -a WholeGenome -v $version
			else
				$java/java -jar $script_path/AddSecondaryAnalysis.jar -p $script_path/AddSecondaryAnalysis.properties -l $lane -f $flowcell -i $index -r $run_num -s Statistics -a WholeGenome -v $version
			fi
			let i=i+1
		done		
	fi
	
	
    if [ $multi_sample != "YES" ]
    then
        echo "Single sample"
		### get the alignment statistics
		if [[ $analysis != "annotation"  && $analysis != "variant" ]] 
		then
			alignment=$input_dir/alignment/$sample	
			cd $alignment
			total_reads=0
			mapped_reads=0
			file=`find . -name '*.flagstat' | sed -e '/\.\//s///g' `
			for i in $file
			do
				re=`cat $i | grep -w 'total'| cut -d ' ' -f1`
				total_reads=`expr $re "+" $total_reads`
			done    
			echo -e "Total Reads" > $numbers/$sample.out
			echo $total_reads >> $numbers/$sample.out
			for i in $file
			do
				ma=` cat $i | grep -w 'mapped' | grep '%)$' | cut -d ' ' -f1`
				mapped_reads=`expr $ma "+" $mapped_reads`
			done    
			echo -e "Mapped Reads (${aligner}) " >> $numbers/$sample.out
			echo $mapped_reads >> $numbers/$sample.out
			percent_dup=0
			echo -e "Percent duplication" >> $numbers/$sample.out
			if [ $markdup == "yes" ]
			then
				## assuming the remove duplicate is done using PICARD
				file=`find . -name '*.dup.metrics' | sed -e '/\.\//s///g' `
				for i in $file
				do
					cat $i | awk '$0 !~ /#/' | cut -f8 | awk '$1 ~ /[^0-9]/' | tail -1 >> $numbers/$sample.dup.out
				done
				percent_dup=`cat $numbers/$sample.dup.out | awk '{sum+=$1; print sum}' | tail -1` 
				rm $numbers/$sample.dup.out 
				echo $percent_dup >> $numbers/$sample.out
			else
				echo $percent_dup >> $numbers/$sample.out		
			fi
		fi
##############################################################		
	    ## statistics after realignment
        if [[ $analysis != "alignment" && $analysis != "annotation" ]]
		then
			realign=$input_dir/realign/$sample
			for chr in $chrs
			do
				rm $realign/chr$chr.cleaned.bam $realign/chr$chr.pileup $realign/chr$chr.cleaned.bam.bai
            done
			cd $realign
            if [ $analysis == "variant" ]
			then
				total_reads=`cat $realign/*.flagstat | grep -w 'total' | cut -f1 | awk '{sum+=$1; print sum}' |tail -1`
				echo -e "Total Reads  " >> $numbers/$sample.out
				echo $total_reads >> $numbers/$sample.out
			fi	
			mapped_reads_realign=0
            mapped_reads_realign=`cat $realign/*.flagstat | grep -w 'mapped' | grep '%)$' | awk '{sum+=$1; print sum}' |tail -1`
            echo -e "Realigned Reads (GATK) " >> $numbers/$sample.out
            echo $mapped_reads_realign >> $numbers/$sample.out
			
			### on target reads
			on_reads=0
			for chr in $chrs
			do
				on=`cat $input_dir/OnTarget/$sample.chr$chr.bam.i.out | head -1`
				on_reads=`expr $on_reads "+" $on`
			done
			echo -e "OnTarget Reads" >> $numbers/$sample.out
			echo $on_reads >> $numbers/$sample.out	
		fi
		if 	[ $analysis != "alignment" ]
		then
            variants=$input_dir/Reports_per_Sample
            ontarget=$input_dir/OnTarget
            ## get the column for the particular sample
            
            if [ $analysis != "annotation" ] 
			then
				## RAW indels and SNVs
				raw_snvs=0
				raw_indels=0
				col=`cat $variants/$sample.variants.raw.vcf | awk '$0 ~ /#/' | tail -1 | awk -v s=$sample -F '\t' '{ for(i=1;i<=NF;i++){ if ($i == s) {print i} } }'`
				a=`cat $variants/$sample.variants.raw.vcf | awk '$0 !~ /#/' | awk -v num=$col '$num !~ /^\./' | awk 'length($4) == 1 && length($5) == 1' | wc -l`
				b=`cat $variants/$sample.variants.raw.vcf | awk '$0 !~ /#/' | awk -v num=$col '$num !~ /^\./' | awk 'length($4) > 1 || length($5) > 1' | wc -l`
				raw_indels=`expr $raw_indels "+" $b`
				raw_snvs=`expr $raw_snvs "+" $a`
				
				
				## Filtered indels and SNVs
				filtered_snvs=0
				filtered_indels=0
				col=`cat $variants/$sample.variants.filter.vcf | awk '$0 ~ /#/' | tail -1 | awk -v s=$sample -F '\t' '{ for(i=1;i<=NF;i++){ if ($i == s) {print i} } }'`
				a=`cat $variants/$sample.variants.filter.vcf | awk '$0 !~ /#/' | awk -v num=$col '$num !~ /^\./' | awk 'length($4) == 1 && length($5) == 1' | grep -c PASS`
				b=`cat $variants/$sample.variants.filter.vcf | awk '$0 !~ /#/' | awk -v num=$col '$num !~ /^\./' | awk 'length($4) > 1 || length($5) > 1' | grep -c PASS`
				filtered_indels=`expr $filtered_indels "+" $b`
				filtered_snvs=`expr $filtered_snvs "+" $a`
			fi
            
            ## Genomic indels and SNVs
			#### SNV in target region anf capture kit
            genomic_snvs=0
            capture_snvs=0
			genomic_indels=0
			capture_indels=0
            for chr in $chrs
            do
                if [ $variant_type == "BOTH" -o $variant_type == "SNV" ]
				then
					s=`cat $ontarget/$sample.chr${chr}.raw.snvs.bed.i.ToMerge | wc -l`
					genomic_snvs=`expr $genomic_snvs "+" $s`
					s_c=`cat $ontarget/$sample.chr${chr}.raw.snvs.bed.i.ToMerge | awk '$(NF-1) == 1' | wc -l`
					capture_snvs=`expr $capture_snvs "+" $s_c`
                fi
				if [ $variant_type == "BOTH" -o $variant_type == "INDEL" ]
				then
					i=`cat $ontarget/$sample.chr${chr}.raw.indels.bed.i.ToMerge | wc -l`
					i_c=`cat $ontarget/$sample.chr${chr}.raw.indels.bed.i.ToMerge | awk '$NF == 1' | wc -l`
					capture_indels=`expr $capture_indels "+" $i_c`
					genomic_indels=`expr $genomic_indels "+" $i`
				fi
			done
            
            if [ $variant_type == "BOTH" -o $variant_type == "SNV" ]
			then
				if [ $analysis != "annotation" ] 
				then
					echo "RAW snvs ($caller)" >> $numbers/$sample.out
					echo $raw_snvs >> $numbers/$sample.out
					echo "FILTERED snvs ($caller)" >> $numbers/$sample.out
					echo $filtered_snvs >> $numbers/$sample.out
				fi
				echo "CODING snvs ($caller)" >> $numbers/$sample.out
				echo $genomic_snvs >> $numbers/$sample.out
			
				if [ $tool == "exome" ]
				then	
					echo "CaptureKit snvs ($caller)" >> $numbers/$sample.out
					echo $capture_snvs >> $numbers/$sample.out	
				fi

				## Annotated SNVs
				sseq_snv=$input_dir/annotation/SSEQ
				touch $sseq_snv/$sample.snv.sseq
				for chr in $chrs
				do
					cat $sseq_snv/$sample.chr${chr}.snv.sseq >> $sseq_snv/$sample.snv.sseq
				done	
				perl $script_path/to.parse.sseq.result.per.sample.pl $sseq_snv/$sample.snv.sseq > $sseq_snv/$sample.snv.sseq.formatted
				
				cat $sseq_snv/$sample.snv.sseq.formatted | awk '$1 ~ "none"' > $sseq_snv/$sample.snv.sseq.formatted.novel
				cat $sseq_snv/$sample.snv.sseq.formatted | awk '$1 !~ "none"' > $sseq_snv/$sample.snv.sseq.formatted.known
				
				# KNOWN variants
				echo -e "Total Known SNVs" >> $numbers/$sample.out
				cat $sseq_snv/$sample.snv.sseq.formatted.known | wc -l >> $numbers/$sample.out
				echo -e "KNOWN Transition To Transversion Ratio" >> $numbers/$sample.out
				perl $script_path/transition.transversion.persample.pl $sseq_snv/$sample.snv.sseq.formatted.known >> $numbers/$sample.out
				echo -e "KNOWN Nonsense" >> $numbers/$sample.out
				cat $sseq_snv/$sample.snv.sseq.formatted.known | awk '$7 ~ "nonsense"' | wc -l >> $numbers/$sample.out
				echo -e "KNOWN Missense" >> $numbers/$sample.out
				cat $sseq_snv/$sample.snv.sseq.formatted.known | awk '$7 ~ "missense"' | wc -l >> $numbers/$sample.out
				echo -e "KNOWN coding-synonymous" >> $numbers/$sample.out
				cat $sseq_snv/$sample.snv.sseq.formatted.known | awk '$7 ~ "coding-synonymous"' | wc -l >> $numbers/$sample.out
				echo -e "KNOWN coding-notMod3" >> $numbers/$sample.out
				cat $sseq_snv/$sample.snv.sseq.formatted.known | awk '$7 ~ "coding-notMod3"' | wc -l >> $numbers/$sample.out
				echo -e "KNOWN splice-3" >> $numbers/$sample.out
				cat $sseq_snv/$sample.snv.sseq.formatted.known | awk '$7 ~ "splice-3"' | wc -l >> $numbers/$sample.out
				echo -e "KNOWN splice-5" >> $numbers/$sample.out
				cat $sseq_snv/$sample.snv.sseq.formatted.known | awk '$7 ~ "splice-5"' | wc -l >> $numbers/$sample.out
				echo -e "KNOWN utr-3" >> $numbers/$sample.out
				cat $sseq_snv/$sample.snv.sseq.formatted.known | awk '$7 ~ "utr-3"' | wc -l >> $numbers/$sample.out
				echo -e "KNOWN utr-5" >> $numbers/$sample.out
				cat $sseq_snv/$sample.snv.sseq.formatted.known | awk '$7 ~ "utr-5"' | wc -l >> $numbers/$sample.out

				# Novel variants
				echo -e "Total Novel SNVs" >> $numbers/$sample.out
				cat $sseq_snv/$sample.snv.sseq.formatted.novel | wc -l >> $numbers/$sample.out
				echo -e "NOVEL Transition To Transversion Ratio" >> $numbers/$sample.out
				perl $script_path/transition.transversion.persample.pl $sseq_snv/$sample.snv.sseq.formatted.novel >> $numbers/$sample.out
				echo -e "NOVEL Nonsense" >> $numbers/$sample.out
				cat $sseq_snv/$sample.snv.sseq.formatted.novel | awk '$7 ~ "nonsense"' | wc -l >> $numbers/$sample.out
				echo -e "NOVEL Missense" >> $numbers/$sample.out
				cat $sseq_snv/$sample.snv.sseq.formatted.novel | awk '$7 ~ "missense"' | wc -l >> $numbers/$sample.out
				echo -e "NOVEL coding-synonymous" >> $numbers/$sample.out
				cat $sseq_snv/$sample.snv.sseq.formatted.novel | awk '$7 ~ "coding-synonymous"' | wc -l >> $numbers/$sample.out
				echo -e "NOVEL coding-notMod3" >> $numbers/$sample.out
				cat $sseq_snv/$sample.snv.sseq.formatted.novel | awk '$7 ~ "coding-notMod3"' | wc -l >> $numbers/$sample.out
				echo -e "NOVEL splice-3" >> $numbers/$sample.out
				cat $sseq_snv/$sample.snv.sseq.formatted.novel | awk '$7 ~ "splice-3"' | wc -l >> $numbers/$sample.out
				echo -e "NOVEL splice-5" >> $numbers/$sample.out
				cat $sseq_snv/$sample.snv.sseq.formatted.novel | awk '$7 ~ "splice-5"' | wc -l >> $numbers/$sample.out
				echo -e "NOVEL utr-3" >> $numbers/$sample.out
				cat $sseq_snv/$sample.snv.sseq.formatted.novel | awk '$7 ~ "utr-3"' | wc -l >> $numbers/$sample.out
				echo -e "NOVEL utr-5" >> $numbers/$sample.out
				cat $sseq_snv/$sample.snv.sseq.formatted.novel | awk '$7 ~ "utr-5"' | wc -l >> $numbers/$sample.out
				
				rm $sseq_snv/$sample.snv.sseq.formatted
				rm $sseq_snv/$sample.snv.sseq.formatted.known
				rm $sseq_snv/$sample.snv.sseq.formatted.novel
				rm $sseq_snv/$sample.snv.sseq
			fi
			
            if [ $variant_type == "BOTH" -o $variant_type == "INDEL" ]
			then
				if [ $analysis != "annotation" ] 
				then
					echo "TOTAL indels ($caller)" >> $numbers/$sample.out
					echo $raw_indels >> $numbers/$sample.out
					echo "FILTERED indels ($caller)" >> $numbers/$sample.out
					echo $filtered_indels >> $numbers/$sample.out
				fi
				echo "CODING indels ($caller)" >> $numbers/$sample.out
				echo $genomic_indels >> $numbers/$sample.out
				if [ $tool == "exome" ]
				then
					echo "Capture indels ($caller)" >> $numbers/$sample.out
					echo $capture_indels >> $numbers/$sample.out
				fi
				## Annotated INDELs
				sseq_indel=$input_dir/annotation/SSEQ
				touch $sseq_indel/$sample.indels.sseq
				for chr in $chrs
				do
					cat $sseq_indel/$sample.chr${chr}.indels.sseq >> $sseq_indel/$sample.indels.sseq
				done	
				perl $script_path/to.parse.sseq.result.indel.per.sample.pl $sseq_indel/$sample.indels.sseq > $sseq_indel/$sample.indels.sseq.formatted
				echo -e "CODING bySSEQ INDELs" >> $numbers/$sample.out
				cat $sseq_indel/$sample.indels.sseq.formatted | awk '$7 ~ "coding"' | wc -l >> $numbers/$sample.out
				echo -e "FRAMESHIFT INDELs" >> $numbers/$sample.out
				cat $sseq_indel/$sample.indels.sseq.formatted | awk '$7 ~ "frameshift"' | wc -l >> $numbers/$sample.out
				echo -e "SPLICE-3 INDELs" >> $numbers/$sample.out
				cat $sseq_indel/$sample.indels.sseq.formatted | awk '$7 ~ "splice-3" ' | wc -l >> $numbers/$sample.out
				echo -e "SPLICE-5 INDELs" >> $numbers/$sample.out
				cat $sseq_indel/$sample.indels.sseq.formatted | awk '$7 ~ "splice-5" ' | wc -l >> $numbers/$sample.out
				echo -e "UTR-3 INDELs" >> $numbers/$sample.out
				cat $sseq_indel/$sample.indels.sseq.formatted | awk '$7 ~ "utr-3" ' | wc -l >> $numbers/$sample.out
				echo -e "UTR-5 INDELs" >> $numbers/$sample.out
				cat $sseq_indel/$sample.indels.sseq.formatted | awk '$7 ~ "utr-5" ' | wc -l >> $numbers/$sample.out
				
				rm $sseq_indel/$sample.indels.sseq
				rm $sseq_indel/$sample.indels.sseq.formatted
			fi
          
			##############################################################		
            if [[ $tool == "whole_genome" && $analysis != "annotation" ]]
			then
				## cnvs
				raw_cnvs=0
				raw_del=0
				raw_dup=0
				genomic_cnvs=0
				genomic_deletions=0
				genomic_duplications=0
				cnv=$input_dir/Reports_per_Sample/SV
				raw_del=`cat $cnv/$sample.cnv.raw.del.vcf | awk '$0 !~ /#/' | wc -l`
				raw_dup=`cat $cnv/$sample.cnv.raw.dup.vcf | awk '$0 !~ /#/' | wc -l`
				raw_cnvs=`expr $raw_del "+" $raw_dup`
				genomic_deletions=`cat $cnv/$sample.cnv.filter.del.vcf | awk '$0 !~ /#/' | wc -l`
				genomic_duplications=`cat $cnv/$sample.cnv.filter.dup.vcf | awk '$0 !~ /#/' | wc -l`
				genomic_cnvs=`expr $genomic_deletions "+" $genomic_duplications`
                rm $cnv/$sample.cnv.raw.del.vcf $cnv/$sample.cnv.raw.dup.vcf $cnv/$sample.cnv.filter.del.vcf $cnv/$sample.cnv.filter.dup.vcf  
				echo "RAW cnvs" >> $numbers/$sample.out
				echo $raw_cnvs >> $numbers/$sample.out
				echo "CODING cnvs" >> $numbers/$sample.out
				echo $genomic_cnvs >> $numbers/$sample.out
				echo "CODING deletions" >> $numbers/$sample.out
				echo $genomic_deletions >> $numbers/$sample.out
				echo "CODING duplications" >> $numbers/$sample.out
				echo $genomic_duplications >> $numbers/$sample.out
            
##############################################################		
				## structural variants 
				struct=$input_dir/Reports_per_Sample/ANNOT
				break=$input_dir/Reports_per_Sample/ANNOT
				
				num_break=0;
				num_crest=0;
				num_sv=0;
				genomic_sv=0;
				ITX=0;
				INV=0;
				DEL=0;
				INS=0;
				CTX=0;
            
				# if [ -d "$struct/crest/$sample" ]
				# then
						# crest=$struct/crest/$sample
				# else
						# echo "CREST folder does not exist"
						# for i in $(seq 1 ${#chrArray[@]})
						# do
								# cat $break/$sample.$chrArray[$i].break.vcf >> $break/$sample.break.vcf
						# done
				num_break=`cat $input_dir/Reports_per_Sample/SV/$sample.break.vcf | awk '$0 !~ /#/' | wc -l`
				num_crest=`cat $input_dir/Reports_per_Sample/SV/$sample.filter.crest.vcf | awk '$0 !~ /#/'|wc -l`
				num_sv=`expr $num_break "+" $num_crest`
				#num_sv=$num_break       
				genomic_sv=`cat $struct/$sample.SV.annotated.txt | wc -l`
				genomic_sv=`expr $genomic_sv "-" 1`
            
				ITX=`cat $struct/$sample.SV.annotated.txt | grep ITX | wc -l`
				INV=`cat $struct/$sample.SV.annotated.txt | grep INV | wc -l`
				DEL=`cat $struct/$sample.SV.annotated.txt | grep DEL | wc -l`
				INS=`cat $struct/$sample.SV.annotated.txt | grep INS | wc -l`
				CTX=`cat $struct/$sample.SV.annotated.txt | grep CTX | wc -l`
				
				echo "TOTAL SVs" >> $numbers/$sample.out
				echo $num_sv >> $numbers/$sample.out
				echo "GENOMIC SVs" >> $numbers/$sample.out
				echo $genomic_sv >> $numbers/$sample.out
				echo "ITX" >> $numbers/$sample.out
				echo $ITX >> $numbers/$sample.out
				echo "INV" >> $numbers/$sample.out
				echo $INV >> $numbers/$sample.out
				echo "DEL" >> $numbers/$sample.out
				echo $DEL >> $numbers/$sample.out
				echo "INS" >> $numbers/$sample.out
				echo $INS >> $numbers/$sample.out
				echo "CTX" >> $numbers/$sample.out
				echo $CTX >> $numbers/$sample.out
			fi
		fi	
	else
	    echo "Multi sample"
        for sample in `cat $sample_info | grep -w "^$group" | cut -d '=' -f2`
        do
            #### alignment stats
            alignment=$input_dir/alignment/$sample	
            cd $alignment
            total_reads=0
            mapped_reads=0
            file=`find . -name '*.flagstat' | sed -e '/\.\//s///g' `
            for i in $file
            do
                re=`cat $i | grep -w 'total'| cut -d ' ' -f1`
                total_reads=`expr $re "+" $total_reads`
            done    
            echo -e "Total Reads" >> $numbers/$sample.out
            echo $total_reads >> $numbers/$sample.out
            for i in $file
            do
                ma=` cat $i | grep -w 'mapped' | grep '%)$' | cut -d ' ' -f1`
                mapped_reads=`expr $ma "+" $mapped_reads`
            done    
            echo -e "Mapped Reads (${aligner}) " >> $numbers/$sample.out
            echo $mapped_reads >> $numbers/$sample.out
            percent_dup=0
            echo -e "Percent duplication" >> $numbers/$sample.out
            if [ $markdup == "yes" ]
            then
                ## assuming the remove duplicate is done using PICARD
                file=`find . -name '*.dup.metrics' | sed -e '/\.\//s///g' `
                for i in $file
                do
                    cat $i | awk '$0 !~ /#/' | cut -f8 | awk '$1 ~ /[^0-9]/' | tail -1 >> $numbers/$sample.dup.out
                done
                percent_dup=`cat $numbers/$sample.dup.out | awk '{sum+=$1; print sum}' | tail -1` 
                rm $numbers/$sample.dup.out 
                echo $percent_dup >> $numbers/$sample.out
            else
                echo $percent_dup >> $numbers/$sample.out		
            fi
            ### on target reads
			on_reads=0
			for chr in $chrs
			do
				on=`cat $input_dir/OnTarget/$sample.chr$chr.bam.i.out | head -1`
				on_reads=`expr $on_reads "+" $on`
			done
			echo -e "OnTarget Reads" >> $numbers/$sample.out
			echo $on_reads >> $numbers/$sample.out	
			
			##### variants
            variants=$input_dir/Reports_per_Sample
            ontarget=$input_dir/OnTarget

            ## RAW indels and SNVs
            raw_snvs=0
            raw_indels=0
            col=`cat $variants/$group.variants.raw.vcf | awk '$0 ~ /#/' | tail -1 | awk -v s=$sample -F '\t' '{ for(i=1;i<=NF;i++){ if ($i == s) {print i} } }'`
            a=`cat $variants/$group.variants.raw.vcf | awk '$0 !~ /#/' | awk -v num=$col '$num !~ /^\./' | awk 'length($4) == 1 && length($5) == 1' | wc -l`
            b=`cat $variants/$group.variants.raw.vcf | awk '$0 !~ /#/' | awk -v num=$col '$num !~ /^\./' | awk 'length($4) > 1 || length($5) > 1' | wc -l`
            raw_indels=`expr $raw_indels "+" $b`
            raw_snvs=`expr $raw_snvs "+" $a`


            ## Filtered indels and SNVs
            filtered_snvs=0
            filtered_indels=0
            col=`cat $variants/$group.variants.filter.vcf | awk '$0 ~ /#/' | tail -1 | awk -v s=$sample -F '\t' '{ for(i=1;i<=NF;i++){ if ($i == s) {print i} } }'`
            a=`cat $variants/$group.variants.filter.vcf | awk '$0 !~ /#/' | awk -v num=$col '$num !~ /^\./' | awk 'length($4) == 1 && length($5) == 1' | grep -c PASS`
            b=`cat $variants/$group.variants.filter.vcf | awk '$0 !~ /#/' | awk -v num=$col '$num !~ /^\./' | awk 'length($4) > 1 || length($5) > 1' | grep -c PASS`
            filtered_indels=`expr $filtered_indels "+" $b`
            filtered_snvs=`expr $filtered_snvs "+" $a`


            ## Genomic indels and SNVs
            #### SNV in target region anf capture kit
            genomic_snvs=0
            genomic_indels=0
            for chr in $chrs
            do
                s=`cat $ontarget/$sample.chr${chr}.raw.snvs.bed.i.ToMerge | wc -l`
                genomic_snvs=`expr $genomic_snvs "+" $s`
                s_c=`cat $ontarget/$sample.chr${chr}.raw.snvs.bed.i.ToMerge | awk '$(NF-1) == 1' | wc -l`
                capture_snvs=`expr $capture_snvs "+" $s_c`
                i=`cat $ontarget/$sample.chr${chr}.raw.indels.bed.i.ToMerge | wc -l`
                i_c=`cat $ontarget/$sample.chr${chr}.raw.indels.bed.i.ToMerge | awk '$NF == 1' | wc -l`
                capture_indels=`expr $capture_indels "+" $i_c`
                genomic_indels=`expr $genomic_indels "+" $i`
            done

            echo "RAW snvs ($caller)" >> $numbers/$sample.out
            echo $raw_snvs >> $numbers/$sample.out
            echo "FILTERED snvs ($caller)" >> $numbers/$sample.out
            echo $filtered_snvs >> $numbers/$sample.out
            echo "CODING snvs ($caller)" >> $numbers/$sample.out
            echo $genomic_snvs >> $numbers/$sample.out
            
            ### annotation 
            sseq_snv=$input_dir/annotation/SSEQ
            touch $sseq_snv/$sample.snv.sseq
            for chr in $chrs
            do
                cat $sseq_snv/$sample.chr${chr}.snv.sseq >> $sseq_snv/$sample.snv.sseq
            done	
            perl $script_path/to.parse.sseq.result.per.sample.pl $sseq_snv/$sample.snv.sseq > $sseq_snv/$sample.snv.sseq.formatted

            cat $sseq_snv/$sample.snv.sseq.formatted | awk '$1 ~ "none"' > $sseq_snv/$sample.snv.sseq.formatted.novel
            cat $sseq_snv/$sample.snv.sseq.formatted | awk '$1 !~ "none"' > $sseq_snv/$sample.snv.sseq.formatted.known

            # KNOWN variants
            echo -e "Total Known SNVs" >> $numbers/$sample.out
            cat $sseq_snv/$sample.snv.sseq.formatted.known | wc -l >> $numbers/$sample.out
            echo -e "KNOWN Transition To Transversion Ratio" >> $numbers/$sample.out
            perl $script_path/transition.transversion.persample.pl $sseq_snv/$sample.snv.sseq.formatted.known >> $numbers/$sample.out
            echo -e "KNOWN Nonsense" >> $numbers/$sample.out
            cat $sseq_snv/$sample.snv.sseq.formatted.known | awk '$7 ~ "nonsense"' | wc -l >> $numbers/$sample.out
            echo -e "KNOWN Missense" >> $numbers/$sample.out
            cat $sseq_snv/$sample.snv.sseq.formatted.known | awk '$7 ~ "missense"' | wc -l >> $numbers/$sample.out
            echo -e "KNOWN coding-synonymous" >> $numbers/$sample.out
            cat $sseq_snv/$sample.snv.sseq.formatted.known | awk '$7 ~ "coding-synonymous"' | wc -l >> $numbers/$sample.out
            echo -e "KNOWN coding-notMod3" >> $numbers/$sample.out
            cat $sseq_snv/$sample.snv.sseq.formatted.known | awk '$7 ~ "coding-notMod3"' | wc -l >> $numbers/$sample.out
            echo -e "KNOWN splice-3" >> $numbers/$sample.out
            cat $sseq_snv/$sample.snv.sseq.formatted.known | awk '$7 ~ "splice-3"' | wc -l >> $numbers/$sample.out
            echo -e "KNOWN splice-5" >> $numbers/$sample.out
            cat $sseq_snv/$sample.snv.sseq.formatted.known | awk '$7 ~ "splice-5"' | wc -l >> $numbers/$sample.out
            echo -e "KNOWN utr-3" >> $numbers/$sample.out
            cat $sseq_snv/$sample.snv.sseq.formatted.known | awk '$7 ~ "utr-3"' | wc -l >> $numbers/$sample.out
            echo -e "KNOWN utr-5" >> $numbers/$sample.out
            cat $sseq_snv/$sample.snv.sseq.formatted.known | awk '$7 ~ "utr-5"' | wc -l >> $numbers/$sample.out

            # Novel variants
            echo -e "Total Novel SNVs" >> $numbers/$sample.out
            cat $sseq_snv/$sample.snv.sseq.formatted.novel | wc -l >> $numbers/$sample.out
            echo -e "NOVEL Transition To Transversion Ratio" >> $numbers/$sample.out
            perl $script_path/transition.transversion.persample.pl $sseq_snv/$sample.snv.sseq.formatted.novel >> $numbers/$sample.out
            echo -e "NOVEL Nonsense" >> $numbers/$sample.out
            cat $sseq_snv/$sample.snv.sseq.formatted.novel | awk '$7 ~ "nonsense"' | wc -l >> $numbers/$sample.out
            echo -e "NOVEL Missense" >> $numbers/$sample.out
            cat $sseq_snv/$sample.snv.sseq.formatted.novel | awk '$7 ~ "missense"' | wc -l >> $numbers/$sample.out
            echo -e "NOVEL coding-synonymous" >> $numbers/$sample.out
            cat $sseq_snv/$sample.snv.sseq.formatted.novel | awk '$7 ~ "coding-synonymous"' | wc -l >> $numbers/$sample.out
            echo -e "NOVEL coding-notMod3" >> $numbers/$sample.out
            cat $sseq_snv/$sample.snv.sseq.formatted.novel | awk '$7 ~ "coding-notMod3"' | wc -l >> $numbers/$sample.out
            echo -e "NOVEL splice-3" >> $numbers/$sample.out
            cat $sseq_snv/$sample.snv.sseq.formatted.novel | awk '$7 ~ "splice-3"' | wc -l >> $numbers/$sample.out
            echo -e "NOVEL splice-5" >> $numbers/$sample.out
            cat $sseq_snv/$sample.snv.sseq.formatted.novel | awk '$7 ~ "splice-5"' | wc -l >> $numbers/$sample.out
            echo -e "NOVEL utr-3" >> $numbers/$sample.out
            cat $sseq_snv/$sample.snv.sseq.formatted.novel | awk '$7 ~ "utr-3"' | wc -l >> $numbers/$sample.out
            echo -e "NOVEL utr-5" >> $numbers/$sample.out
            cat $sseq_snv/$sample.snv.sseq.formatted.novel | awk '$7 ~ "utr-5"' | wc -l >> $numbers/$sample.out

            echo "TOTAL indels ($caller)" >> $numbers/$sample.out
            echo $raw_indels >> $numbers/$sample.out
            echo "FILTERED indels ($caller)" >> $numbers/$sample.out
            echo $filtered_indels >> $numbers/$sample.out
            echo "CODING indels ($caller)" >> $numbers/$sample.out
            echo $genomic_indels >> $numbers/$sample.out
            sseq_indel=$input_dir/annotation/SSEQ
            touch $sseq_indel/$sample.indels.sseq
            for chr in $chrs
            do
                cat $sseq_indel/$sample.chr${chr}.indels.sseq >> $sseq_indel/$sample.indels.sseq
            done	
            perl $script_path/to.parse.sseq.result.indel.per.sample.pl $sseq_indel/$sample.indels.sseq > $sseq_indel/$sample.indels.sseq.formatted
            echo -e "CODING bySSEQ INDELs" >> $numbers/$sample.out
            cat $sseq_indel/$sample.indels.sseq.formatted | awk '$7 ~ "coding"' | wc -l >> $numbers/$sample.out
            echo -e "FRAMESHIFT INDELs" >> $numbers/$sample.out
            cat $sseq_indel/$sample.indels.sseq.formatted | awk '$7 ~ "frameshift"' | wc -l >> $numbers/$sample.out
            echo -e "SPLICE-3 INDELs" >> $numbers/$sample.out
            cat $sseq_indel/$sample.indels.sseq.formatted | awk '$7 ~ "splice-3" ' | wc -l >> $numbers/$sample.out
            echo -e "SPLICE-5 INDELs" >> $numbers/$sample.out
            cat $sseq_indel/$sample.indels.sseq.formatted | awk '$7 ~ "splice-5" ' | wc -l >> $numbers/$sample.out
            echo -e "UTR-3 INDELs" >> $numbers/$sample.out
            cat $sseq_indel/$sample.indels.sseq.formatted | awk '$7 ~ "utr-3" ' | wc -l >> $numbers/$sample.out
            echo -e "UTR-5 INDELs" >> $numbers/$sample.out
            cat $sseq_indel/$sample.indels.sseq.formatted | awk '$7 ~ "utr-5" ' | wc -l >> $numbers/$sample.out

            rm $sseq_snv/$sample.snv.sseq.formatted
            rm $sseq_snv/$sample.snv.sseq.formatted.known
            rm $sseq_snv/$sample.snv.sseq.formatted.novel
            rm $sseq_indel/$sample.indels.sseq
            rm $sseq_indel/$sample.indels.sseq.formatted
            rm $sseq_snv/$sample.snv.sseq
        done
        
		### somatic calls
        sample_info=$( cat $run_info | grep -w '^SAMPLE_INFO' | cut -d '=' -f2)
        samples=$( cat $sample_info | grep -w "^$group" | cut -d '=' -f2 )
        let num_tumor=`echo $samples|tr " " "\n"|wc -l`-1
        tumor_list=`echo $samples | tr " " "\n" | tail -$num_tumor`
        for tumor in $tumor_list    
        do
            realignment=$input_dir/realign/$group	
            total_reads=0
            mapped_reads=0
            for chr in $chrs
            do
                re=`cat $realignment/chr$chr.flagstat | grep -w 'total'| cut -d ' ' -f1`
                total_reads=`expr $re "+" $total_reads`
            done    
            echo -e "Combined Total Reads After alignment" > $numbers/$group.$tumor.out
            echo $total_reads >> $numbers/$group.$tumor.out
            for chr in $chrs
            do
                ma=` cat $realignment/chr$chr.flagstat | grep -w 'mapped' | grep '%)$' | cut -d ' ' -f1`
                mapped_reads=`expr $ma "+" $mapped_reads`
            done    
            echo -e "Mapped Reads (GATK) " >> $numbers/$group.$tumor.out
            echo $mapped_reads >> $numbers/$group.$tumor.out
            ##### variants
            variants=$input_dir/Reports_per_Sample
            ontarget=$input_dir/OnTarget

            ## RAW indels and SNVs
            raw_snvs=0
            raw_indels=0
            col=`cat $variants/$group.somatic.variants.raw.vcf | awk '$0 ~ /#/' | tail -1 | awk -v s=$sample -F '\t' '{ for(i=1;i<=NF;i++){ if ($i == s) {print i} } }'`
            a=`cat $variants/$group.somatic.variants.raw.vcf | awk '$0 !~ /#/' | awk -v num=$col '$num !~ /^\./' | awk 'length($4) == 1 && length($5) == 1' | wc -l`
            b=`cat $variants/$group.somatic.variants.raw.vcf | awk '$0 !~ /#/' | awk -v num=$col '$num !~ /^\./' | awk 'length($4) > 1 || length($5) > 1' | wc -l`
            raw_indels=`expr $raw_indels "+" $b`
            raw_snvs=`expr $raw_snvs "+" $a`


            ## Filtered indels and SNVs
            filtered_snvs=0
            filtered_indels=0
            col=`cat $variants/$group.somatic.variants.filter.vcf | awk '$0 ~ /#/' | tail -1 | awk -v s=$sample -F '\t' '{ for(i=1;i<=NF;i++){ if ($i == s) {print i} } }'`
            a=`cat $variants/$group.somatic.variants.filter.vcf | awk '$0 !~ /#/' | awk -v num=$col '$num !~ /^\./' | awk 'length($4) == 1 && length($5) == 1' | wc -l`
            b=`cat $variants/$group.somatic.variants.filter.vcf | awk '$0 !~ /#/' | awk -v num=$col '$num !~ /^\./' | awk 'length($4) > 1 || length($5) > 1' | wc -l `
            filtered_indels=`expr $filtered_indels "+" $b`
            filtered_snvs=`expr $filtered_snvs "+" $a`


            ## Genomic indels and SNVs
            #### SNV in target region anf capture kit
            genomic_snvs=0
            capture_snvs=0
            genomic_indels=0
            capture_indels=0
            for chr in $chrs
            do
                s=`cat $ontarget/$group.$tumor.chr${chr}.raw.snvs.bed.i.ToMerge | wc -l`
                genomic_snvs=`expr $genomic_snvs "+" $s`
                s_c=`cat $ontarget/$group.$tumor.chr${chr}.raw.snvs.bed.i.ToMerge | awk '$(NF-1) == 1' | wc -l`
                capture_snvs=`expr $capture_snvs "+" $s_c`
                i=`cat $ontarget/$group.$tumor.chr${chr}.raw.indels.bed.i.ToMerge | wc -l`
                i_c=`cat $ontarget/$group.$tumor.chr${chr}.raw.indels.bed.i.ToMerge | awk '$NF == 1' | wc -l`
                capture_indels=`expr $capture_indels "+" $i_c`
                genomic_indels=`expr $genomic_indels "+" $i`
            done

            echo "RAW snvs ($caller)" >> $numbers/$group.$tumor.out
            echo $raw_snvs >> $numbers/$group.$tumor.out
            echo "FILTERED snvs ($caller)" >> $numbers/$group.$tumor.out
            echo $filtered_snvs >> $numbers/$group.$tumor.out
            echo "CODING snvs ($caller)" >> $numbers/$group.$tumor.out
            echo $genomic_snvs >> $numbers/$group.$tumor.out

            ### annotation 
            sseq_snv=$input_dir/annotation/SSEQ
            touch $sseq_snv/$group.$tumor.snv.sseq
            for chr in $chrs
            do
                cat $sseq_snv/$group.$tumor.chr${chr}.snv.sseq >> $sseq_snv/$group.$tumor.snv.sseq
            done	
            perl $script_path/to.parse.sseq.result.per.sample.pl $sseq_snv/$group.$tumor.snv.sseq > $sseq_snv/$group.$tumor.snv.sseq.formatted

            cat $sseq_snv/$group.$tumor.snv.sseq.formatted | awk '$1 ~ "none"' > $sseq_snv/$group.$tumor.snv.sseq.formatted.novel
            cat $sseq_snv/$group.$tumor.snv.sseq.formatted | awk '$1 !~ "none"' > $sseq_snv/$group.$tumor.snv.sseq.formatted.known

            # KNOWN variants
            echo -e "Total Known SNVs" >> $numbers/$group.$tumor.out
            cat $sseq_snv/$group.$tumor.snv.sseq.formatted.known | wc -l >> $numbers/$group.$tumor.out
            echo -e "KNOWN Transition To Transversion Ratio" >> $numbers/$group.$tumor.out
            perl $script_path/transition.transversion.persample.pl $sseq_snv/$group.$tumor.snv.sseq.formatted.known >> $numbers/$group.$tumor.out
            echo -e "KNOWN Nonsense" >> $numbers/$group.$tumor.out
            cat $sseq_snv/$group.$tumor.snv.sseq.formatted.known | awk '$7 ~ "nonsense"' | wc -l >> $numbers/$group.$tumor.out
            echo -e "KNOWN Missense" >> $numbers/$group.$tumor.out
            cat $sseq_snv/$group.$tumor.snv.sseq.formatted.known | awk '$7 ~ "missense"' | wc -l >> $numbers/$group.$tumor.out
            echo -e "KNOWN coding-synonymous" >> $numbers/$group.$tumor.out
            cat $sseq_snv/$group.$tumor.snv.sseq.formatted.known | awk '$7 ~ "coding-synonymous"' | wc -l >> $numbers/$group.$tumor.out
            echo -e "KNOWN coding-notMod3" >> $numbers/$group.$tumor.out
            cat $sseq_snv/$group.$tumor.snv.sseq.formatted.known | awk '$7 ~ "coding-notMod3"' | wc -l >> $numbers/$group.$tumor.out
            echo -e "KNOWN splice-3" >> $numbers/$group.$tumor.out
            cat $sseq_snv/$group.$tumor.snv.sseq.formatted.known | awk '$7 ~ "splice-3"' | wc -l >> $numbers/$group.$tumor.out
            echo -e "KNOWN splice-5" >> $numbers/$group.$tumor.out
            cat $sseq_snv/$group.$tumor.snv.sseq.formatted.known | awk '$7 ~ "splice-5"' | wc -l >> $numbers/$group.$tumor.out
            echo -e "KNOWN utr-3" >> $numbers/$group.$tumor.out
            cat $sseq_snv/$group.$tumor.snv.sseq.formatted.known | awk '$7 ~ "utr-3"' | wc -l >> $numbers/$group.$tumor.out
            echo -e "KNOWN utr-5" >> $numbers/$group.$tumor.out
            cat $sseq_snv/$group.$tumor.snv.sseq.formatted.known | awk '$7 ~ "utr-5"' | wc -l >> $numbers/$group.$tumor.out

            # Novel variants
            echo -e "Total Novel SNVs" >> $numbers/$group.$tumor.out
            cat $sseq_snv/$group.$tumor.snv.sseq.formatted.novel | wc -l >> $numbers/$group.$tumor.out
            echo -e "NOVEL Transition To Transversion Ratio" >> $numbers/$group.$tumor.out
            perl $script_path/transition.transversion.persample.pl $sseq_snv/$group.$tumor.snv.sseq.formatted.novel >> $numbers/$group.$tumor.out
            echo -e "NOVEL Nonsense" >> $numbers/$group.$tumor.out
            cat $sseq_snv/$group.$tumor.snv.sseq.formatted.novel | awk '$7 ~ "nonsense"' | wc -l >> $numbers/$group.$tumor.out
            echo -e "NOVEL Missense" >> $numbers/$group.$tumor.out
            cat $sseq_snv/$group.$tumor.snv.sseq.formatted.novel | awk '$7 ~ "missense"' | wc -l >> $numbers/$group.$tumor.out
            echo -e "NOVEL coding-synonymous" >> $numbers/$group.$tumor.out
            cat $sseq_snv/$group.$tumor.snv.sseq.formatted.novel | awk '$7 ~ "coding-synonymous"' | wc -l >> $numbers/$group.$tumor.out
            echo -e "NOVEL coding-notMod3" >> $numbers/$group.$tumor.out
            cat $sseq_snv/$group.$tumor.snv.sseq.formatted.novel | awk '$7 ~ "coding-notMod3"' | wc -l >> $numbers/$group.$tumor.out
            echo -e "NOVEL splice-3" >> $numbers/$group.$tumor.out
            cat $sseq_snv/$group.$tumor.snv.sseq.formatted.novel | awk '$7 ~ "splice-3"' | wc -l >> $numbers/$group.$tumor.out
            echo -e "NOVEL splice-5" >> $numbers/$group.$tumor.out
            cat $sseq_snv/$group.$tumor.snv.sseq.formatted.novel | awk '$7 ~ "splice-5"' | wc -l >> $numbers/$group.$tumor.out
            echo -e "NOVEL utr-3" >> $numbers/$group.$tumor.out
            cat $sseq_snv/$group.$tumor.snv.sseq.formatted.novel | awk '$7 ~ "utr-3"' | wc -l >> $numbers/$group.$tumor.out
            echo -e "NOVEL utr-5" >> $numbers/$group.$tumor.out
            cat $sseq_snv/$group.$tumor.snv.sseq.formatted.novel | awk '$7 ~ "utr-5"' | wc -l >> $numbers/$group.$tumor.out

            echo "TOTAL indels ($caller)" >> $numbers/$group.$tumor.out
            echo $raw_indels >> $numbers/$group.$tumor.out
            echo "FILTERED indels ($caller)" >> $numbers/$group.$tumor.out
            echo $filtered_indels >> $numbers/$group.$tumor.out
            echo "CODING indels ($caller)" >> $numbers/$group.$tumor.out
            echo $genomic_indels >> $numbers/$group.$tumor.out
            sseq_indel=$input_dir/annotation/SSEQ
            touch $sseq_indel/$group.$tumor.indels.sseq
            for chr in $chrs
            do
                cat $sseq_indel/$group.$tumor.chr${chr}.indels.sseq >> $sseq_indel/$group.$tumor.indels.sseq
            done	
            perl $script_path/to.parse.sseq.result.indel.per.sample.pl $sseq_indel/$group.$tumor.indels.sseq > $sseq_indel/$group.$tumor.indels.sseq.formatted
            echo -e "CODING bySSEQ INDELs" >> $numbers/$group.$tumor.out
            cat $sseq_indel/$group.$tumor.indels.sseq.formatted | awk '$7 ~ "coding"' | wc -l >> $numbers/$group.$tumor.out
            echo -e "FRAMESHIFT INDELs" >> $numbers/$group.$tumor.out
            cat $sseq_indel/$group.$tumor.indels.sseq.formatted | awk '$7 ~ "frameshift"' | wc -l >> $numbers/$group.$tumor.out
            echo -e "SPLICE-3 INDELs" >> $numbers/$group.$tumor.out
            cat $sseq_indel/$group.$tumore.indels.sseq.formatted | awk '$7 ~ "splice-3" ' | wc -l >>$numbers/$group.$tumor.out
            echo -e "SPLICE-5 INDELs" >> $numbers/$group.$tumor.out
            cat $sseq_indel/$group.$tumor.indels.sseq.formatted | awk '$7 ~ "splice-5" ' | wc -l >> $numbers/$group.$tumor.out
            echo -e "UTR-3 INDELs" >> $numbers/$group.$tumor.out
            cat $sseq_indel/$group.$tumor.indels.sseq.formatted | awk '$7 ~ "utr-3" ' | wc -l >> $numbers/$group.$tumor.out
            echo -e "UTR-5 INDELs" >> $numbers/$group.$tumor.out
            cat $sseq_indel/$group.$tumor.indels.sseq.formatted | awk '$7 ~ "utr-5" ' | wc -l >> $numbers/$group.$tumor.out

            rm $sseq_snv/$group.$tumor.snv.sseq.formatted
            rm $sseq_snv/$group.$tumor.snv.sseq.formatted.known
            rm $sseq_snv/$group.$tumor.snv.sseq.formatted.novel
            rm $sseq_indel/$group.$tumor.indels.sseq
            rm $sseq_indel/$group.$tumor.indels.sseq.formatted
            rm $sseq_snv/$group.$tumor.snv.sseq
            ###   cnvs by segseq
            raw_cnvs=0
            raw_del=0
            raw_dup=0
            genomic_cnvs=0
            genomic_deletions=0
            genomic_duplications=0
            cnv=$input_dir/Reports_per_Sample/SV
            raw_del=`cat $cnv/$group.$tumor.cnv.raw.del.vcf | awk '$0 !~ /#/' | wc -l`
            raw_dup=`cat $cnv/$group.$tumor.cnv.raw.dup.vcf | awk '$0 !~ /#/' | wc -l`
            raw_cnvs=`expr $raw_del "+" $raw_dup`
            genomic_deletions=`cat $cnv/$group.$tumor.cnv.filter.del.vcf | awk '$0 !~ /#/' | wc -l`
            genomic_duplications=`cat $cnv/$group.$tumor.cnv.filter.dup.vcf | awk '$0 !~ /#/' | wc -l`
            genomic_cnvs=`expr $genomic_deletions "+" $genomic_duplications`
            rm $cnv/$group.$tumor.cnv.raw.del.vcf $cnv/$group.$tumor.cnv.raw.dup.vcf $cnv/$group.$tumor.cnv.filter.del.vcf $cnv/$group.$tumor.cnv.filter.dup.vcf
            echo "RAW cnvs" >> $numbers/$group.$tumor.out
            echo $raw_cnvs >> $numbers/$group.$tumor.out
            echo "CODING cnvs" >> $numbers/$group.$tumor.out
            echo $genomic_cnvs >> $numbers/$group.$tumor.out
            echo "CODING deletions" >> $numbers/$group.$tumor.out
            echo $genomic_deletions >> $numbers/$group.$tumor.out
            echo "CODING duplications" >> $numbers/$group.$tumor.out
            echo $genomic_duplications >> $numbers/$group.$tumor.out
            #### SV by crest and break dancer
			struct=$input_dir/Reports_per_Sample/ANNOT
			break=$input_dir/Reports_per_Sample/ANNOT
				
			num_break=0;
			num_crest=0;
			num_sv=0;
			genomic_sv=0;
			ITX=0;
			INV=0;
			DEL=0;
			INS=0;
			CTX=0;
            num_break=`cat $input_dir/Reports_per_Sample/SV/$group.$tumor.somatic.break.vcf | awk '$0 !~ /#/' | wc -l`
			num_crest=`cat $input_dir/Reports_per_Sample/SV/$group.$tumor.somatic.filter.crest.vcf | awk '$0 !~ /#/'|wc -l`
			num_sv=`expr $num_break "+" $num_crest`
			genomic_sv=`cat $struct/$group.$tumor.SV.annotated.txt | wc -l`
			genomic_sv=`expr $genomic_sv "-" 1`
		
			ITX=`cat $struct/$group.$tumor.SV.annotated.txt | grep ITX | wc -l`
			INV=`cat $struct/$group.$tumor.SV.annotated.txt | grep INV | wc -l`
			DEL=`cat $struct/$group.$tumor.SV.annotated.txt | grep DEL | wc -l`
			INS=`cat $struct/$group.$tumor.SV.annotated.txt | grep INS | wc -l`
			CTX=`cat $struct/$group.$tumor.SV.annotated.txt | grep CTX | wc -l`
			
			echo "TOTAL SVs" >> $numbers/$group.$tumor.out
			echo $num_sv >> $numbers/$group.$tumor.out
			echo "GENOMIC SVs" >> $numbers/$group.$tumor.out
			echo $genomic_sv >> $numbers/$group.$tumor.out
			echo "ITX" >>$numbers/$group.$tumor.out
			echo $ITX >> $numbers/$group.$tumor.out
			echo "INV" >> $numbers/$group.$tumor.out
			echo $INV >> $numbers/$group.$tumor.out
			echo "DEL" >> $numbers/$group.$tumor.out
			echo $DEL >> $numbers/$group.$tumor.out
			echo "INS" >> $numbers/$group.$tumor.out
			echo $INS >> $numbers/$group.$tumor.out
			echo "CTX" >> $numbers/$group.$tumor.out
			echo $CTX >> $numbers/$group.$tumor.out
        done
	fi 
	### update dash board
	if [ $analysis == "mayo" -o $analysis == "realign-mayo" ]
	then
		pos=$( cat $run_info | grep -w '^SAMPLENAMES' | cut -d '=' -f2 | tr ":" "\n" | grep -n $sample | cut -d ":" -f1)
		lanes=$( cat $run_info | grep -w '^LANEINDEX' | cut -d '=' -f2 | tr ":" "\n" | head -n $pos | tail -n 1 | tr "," " ")
		i=1
		for lane in $lanes
		do
			index=$( cat $run_info | grep -w '^LABINDEXES' | cut -d '=' -f2 | tr ":" "\n" | head -n $pos | tail -n 1 | tr "," "\n" | head -n $i | tail -n 1)
			if [ $index == "-" ]
			then
				$java/java -jar $script_path/AddSecondaryAnalysis.jar -p $script_path/AddSecondaryAnalysis.properties -l $lane -c -f $flowcell -r $run_num -s Statistics -a WholeGenome -v $version
			else
				$java/java -jar $script_path/AddSecondaryAnalysis.jar -p $script_path/AddSecondaryAnalysis.properties -l $lane -c -f $flowcell -i $index -r $run_num -s Statistics -a WholeGenome -v $version
			fi
			let i=i+1
		done		
	fi
	echo `date`
fi   
