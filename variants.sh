#!/bin/bash

if [ $# -le 4 ]
then
	echo -e "script to run variant calling on set of BAM files\nUsage: ./variants.sh </path/to/input dir><samples ':' sep[normal:tumor1:tumor2:tumorN]> </path/to/output directory> <1 or 0 chopped or not ><path/to/run info><SGE_TASK_ID(optional)>\nNOTE: first sample is considered as normal and others are considered as tumor[Assumption]\n";
else
	set -x
	echo `date`
	input=$1
	samples=$2
	output=$3
	chopped=$4
	run_info=$5
	if [ $6 ]
	then
		SGE_TASK_ID=$6
	fi	
	tool_info=$( cat $run_info | grep -w '^TOOL_INFO' | cut -d '=' -f2)
	tool=$( cat $run_info | grep -w '^TYPE' | cut -d '=' -f2|tr "[A-Z]" "[a-z]")
	tabix=$( cat $tool_info | grep -w '^TABIX' | cut -d '=' -f2)
	perllib=$( cat $tool_info | grep -w '^PERLLIB_VCF' | cut -d '=' -f2)
	samtools=$( cat $tool_info | grep -w '^SAMTOOLS' | cut -d '=' -f2)
	ref=$( cat $tool_info | grep -w '^REF_GENOME' | cut -d '=' -f2)
	dbSNP=$( cat $tool_info | grep -w '^dbSNP_REF' | cut -d '=' -f2)
	Kgenome=$( cat $tool_info | grep -w '^KGENOME_REF' | cut -d '=' -f2)
	java=$( cat $tool_info | grep -w '^JAVA' | cut -d '=' -f2)
	picard=$( cat $tool_info | grep -w '^PICARD' | cut -d '=' -f2 )
	script_path=$( cat $tool_info | grep -w '^WORKFLOW_PATH' | cut -d '=' -f2 )
	analysis=$( cat $run_info | grep -w '^ANALYSIS' | cut -d '=' -f2 |tr "[A-Z]" "[a-z]")
	sampleNames=$( echo $samples | tr ":" "\n" )
	chr=$( cat $run_info | grep -w '^CHRINDEX' | cut -d '=' -f2 | tr ":" "\n" | head -n $SGE_TASK_ID | tail -n 1 )
	all_sites=$( cat $tool_info | grep -w '^EMIT_ALL_SITES' | cut -d '=' -f2 | tr "[a-z]" "[A-Z]" )
	depth_filter=$( cat $tool_info | grep -w '^DEPTH_FILTER' | cut -d '=' -f2)
	only_ontarget=$( cat $tool_info | grep -w '^TARGETTED' | cut -d '=' -f2 | tr "[a-z]" "[A-Z]" )
	TargetKit=$( cat $tool_info | grep -w '^ONTARGET' | cut -d '=' -f2 )
	SNV_caller=$( cat $run_info | grep -w '^SNV_CALLER' | cut -d '=' -f2)
	somatic_caller=$( cat $run_info | grep -w '^SOMATIC_CALLER' | cut -d '=' -f2)
	ped=$( cat $tool_info | grep -w '^PEDIGREE' | cut -d '=' -f2 )
	bedtools=$( cat $tool_info | grep -w '^BEDTOOLS' | cut -d '=' -f2 )
	blat=$( cat $tool_info | grep -w '^BLAT' | cut -d '=' -f2 )
	blat_ref=$( cat $tool_info | grep -w '^BLAT_REF' | cut -d '=' -f2 )
	blat_params=$( cat $tool_info | grep -w '^BLAT_params' | cut -d '=' -f2 )
	perlblat=$( cat $tool_info | grep -w '^PERLLIB_BLAT' | cut -d '=' -f2 )
	export PERL5LIB=$PERL5LIB:$perllib:$perlblat
	export PATH=$java:$tabix/:$PATH
	
	#### check and validate the bam file and let user to proceed after validation
	bam=chr${chr}.cleaned.bam
	$samtools/samtools view -H $input/$bam 1>$input/$bam.$chr.header 2> $input/$bam.$chr.fix.log
	if [ `cat $input/$bam.$chr.fix.log | wc -l` -gt 0 ]
	then
		$script_path/email.sh $input/$bam "bam is truncated or corrupt" realign_recal.sh $run_info
		$script_path/wait.sh $input/$bam.$chr.fix.log
	else
		rm $input/$bam.$chr.fix.log
	fi		
	rm $input/$bam.$chr.header	
	if [ `echo $samples | tr ":" "\n" | wc -l` -gt 1 ]
	then
		$script_path/filesize.sh VariantCalling multi_sample $input $bam $run_info
	else
		$script_path/filesize.sh VariantCalling $samples $input $bam $run_info
	fi
	
	## update dashborad
	if [ $SGE_TASK_ID == 1 ]
	then
		for i in `echo $samples | tr ":" " "`
		do
			$script_path/dashboard.sh $i $run_info VariantCalling started
		done
	fi

	i=1
	for sample in $sampleNames
	do
		sampleArray[$i]=$sample
		let i=i+1
	done
	num_samples=`echo $samples | tr ":" " "`

	mkdir -p $output/temp
	## if multiple samples then split using read group and Validate BAM other wise just check and validate BAM
	if [ ${#sampleArray[@]} -gt 1 ]
	then
		## bams are splitted using read group information
		for i in $(seq 1 ${#sampleArray[@]})
		do
			sample=${sampleArray[$i]}
			### removing the header for the extra samples from the BAM
			sam=`echo $samples | tr ":" "\n"| grep -w -v "$sample" | tr "\n" " "`
			gr=""
			for s in $sam
			do
				a="ID:$s|";
				gr="$gr$a"
			done
			gr=`echo $gr |  sed "s/|$//"`
			$samtools/samtools view -b -r $sample $input/$bam > $output/$sample.chr$chr.rg.bam
			$samtools/samtools view -H $output/$sample.chr$chr.rg.bam | grep -w -E -v "$gr" | $samtools/samtools reheader - $output/$sample.chr$chr.rg.bam > $output/$sample.chr$chr.rg.re.bam
			mv $output/$sample.chr$chr.rg.re.bam $output/$sample.chr$chr.rg.bam

			if [ ! -s $output/$sample.chr$chr.rg.bam ]
			then
				$script_path/errorlog.sh $output/$sample.chr$chr.rg.bam variants.sh ERROR "failed to create"
				exit 1;
			fi
			$script_path/samplecheckBAM.sh $output $sample.chr$chr.rg.bam $output $run_info $sample $chopped $chr
		done
	else
		sample=${sampleArray[1]}
		$script_path/samplecheckBAM.sh $input $bam $output $run_info $sample $chopped $chr
	fi

	inputfiles=""
	if [ $tool == "exome" ]
	then
		cat $TargetKit | grep -w chr$chr > $output/chr$chr.target.bed
	fi
	if [ ${#sampleArray[@]} == 1 ]
	then
		sample=${sampleArray[1]}
		if [ $SNV_caller == "GATK" ]
		then
			## call variants using gatk UnifiedGenotyper module
			if [[ $all_sites == "YES"  && $tool == "exome" ]]
			then
				len=`cat $output/chr$chr.target.bed |wc -l`
				param="-L chr${chr}"
				bam="-I $output/$sample.chr${chr}-sorted.bam"
				$script_path/unifiedgenotyper.sh "$bam" $output/$sample.variants.chr${chr}.raw.all.vcf BOTH "$param" EMIT_ALL_SITES $run_info
				## add phase by transmission if pedigree information provided.
				if [ $ped != "NA" ]
				then
					$script_path/phaseByTransmission.sh $output/$sample.variants.chr${chr}.raw.all.vcf $output/$sample.variants.chr${chr}.raw.all.pbt.vcf $run_info
				fi
				### filter this file to keep only the variants calls
				cat $output/$sample.variants.chr${chr}.raw.all.vcf | awk '$5 != "." || $0 ~ /^#/' | grep -v "\./\."  | grep -v "0\/0" > $output/$sample.variants.chr${chr}.raw.vcf
				sed '/^$/d' $output/$sample.variants.chr${chr}.raw.vcf > $output/$sample.variants.chr${chr}.raw.vcf.temp
				mv $output/$sample.variants.chr${chr}.raw.vcf.temp $output/$sample.variants.chr${chr}.raw.vcf
				### prepare the file for backfilling
				cat $output/$sample.variants.chr${chr}.raw.all.vcf | grep -v "\./\." > $output/$sample.variants.chr${chr}.raw.all.vcf.temp
				mv $output/$sample.variants.chr${chr}.raw.all.vcf.temp $output/$sample.variants.chr${chr}.raw.all.vcf
				sed '/^$/d' $output/$sample.variants.chr${chr}.raw.all.vcf > $output/$sample.variants.chr${chr}.raw.all.vcf.temp
				mv $output/$sample.variants.chr${chr}.raw.all.vcf.temp $output/$sample.variants.chr${chr}.raw.all.vcf
				if [[ $len -gt 0 && $only_ontarget == "YES" ]]
				then
					$bedtools/intersectBed -a $output/$sample.variants.chr${chr}.raw.all.vcf -b $output/chr$chr.target.bed -wa -header > $output/$sample.variants.chr${chr}.raw.all.vcf.temp
					mv $output/$sample.variants.chr${chr}.raw.all.vcf.temp $output/$sample.variants.chr${chr}.raw.all.vcf
				fi
				$tabix/bgzip $output/$sample.variants.chr${chr}.raw.all.vcf
				mv $output/$sample.variants.chr${chr}.raw.all.vcf.multi.vcf $output/$sample.variants.chr${chr}.raw.multi.vcf
			else
				param="-L chr${chr}"
				bam="-I $output/$sample.chr${chr}-sorted.bam"
				$script_path/unifiedgenotyper.sh "$bam" $output/$sample.variants.chr${chr}.raw.vcf BOTH "$param" EMIT_VARIANTS_ONLY $run_info
				## add phase by transmission if pedigree information provided.
				if [ $ped != "NA" ]
				then
					$script_path/phaseByTransmission.sh $output/$sample.variants.chr${chr}.raw.all.vcf $output/$sample.variants.chr${chr}.raw.all.pbt.vcf $run_info
				fi
				mv $output/$sample.variants.chr${chr}.raw.vcf.multi.vcf $output/$sample.variants.chr${chr}.raw.multi.vcf
				$script_path/annotate_vcf.sh $output/$sample.variants.chr${chr}.raw.multi.vcf $run_info "$bam" 
				$script_path/annotate_vcf.sh $output/$sample.variants.chr${chr}.raw.vcf $run_info "$bam" 
			fi
		elif [ $SNV_caller == "SNVMIX" ]
		then
			### call indels using GATK
			if [[ $all_sites == "YES"  && $tool == "exome" ]]
			then
				len=`cat $output/chr$chr.target.bed |wc -l`
				param="-L chr${chr}"
				## call indels
				bam="-I $output/$sample.chr${chr}-sorted.bam"
				$script_path/unifiedgenotyper.sh "$bam" $output/$sample.variants.chr${chr}.raw.indel.all.vcf INDEL "$param" EMIT_ALL_SITES $run_info
				### call snvs using SNVmix
				$script_path/snvmix2.sh $sample "$bam" $output/$sample.variants.chr${chr}.raw.snv.all.vcf all "$param" $run_info
				### annoatte vcf
				$script_path/annotate_vcf.sh $output/$sample.variants.chr${chr}.raw.snv.all.vcf $run_info "$bam"
				### merge snvs and indels to give on vcf
				in="$output/$sample.variants.chr${chr}.raw.snv.all.vcf $output/$sample.variants.chr${chr}.raw.indel.all.vcf "
				$script_path/concatvcf.sh "$in" $output/$sample.variants.chr${chr}.raw.all.vcf $run_info yes
				in="$output/$sample.variants.chr${chr}.raw.snv.all.vcf.multi.vcf $output/$sample.variants.chr${chr}.raw.indel.all.vcf.multi.vcf"
				$script_path/concatvcf.sh "$in" $output/$sample.variants.chr${chr}.raw.multi.vcf $run_info yes
				cat $output/$sample.variants.chr${chr}.raw.all.vcf | awk '$5 != "." || $0 ~ /^#/' | grep -v "\./\."  | grep -v "0\/0" > $output/$sample.variants.chr${chr}.raw.vcf
				if [[ $len -gt 0 && $only_ontarget == "YES" ]]
				then
					$bedtools/intersectBed -a $output/$sample.variants.chr${chr}.raw.all.vcf -b $output/chr$chr.target.bed -wa -header > $output/$sample.variants.chr${chr}.raw.all.vcf.temp
					mv $output/$sample.variants.chr${chr}.raw.all.vcf.temp $output/$sample.variants.chr${chr}.raw.all.vcf
				fi
				$tabix/bgzip $output/$sample.variants.chr${chr}.raw.all.vcf
			else
				## call indeles using GATK
				param="-L chr${chr}"
				bam="-I $output/$sample.chr${chr}-sorted.bam"
				$script_path/unifiedgenotyper.sh "$bam" $output/$sample.variants.chr${chr}.raw.indel.vcf INDEL "$param" EMIT_VARIANTS_ONLY $run_info &
				### call snvs using snvmix
				$script_path/snvmix2.sh $sample "$bam" $output/$sample.variants.chr${chr}.raw.snv.vcf target "$param" $run_info &
				while [[ ! -s $output/$sample.variants.chr${chr}.raw.indel.vcf || ! -s $output/$sample.variants.chr${chr}.raw.indel.vcf.multi.vcf || ! -s $output/$sample.variants.chr${chr}.raw.snv.vcf || ! -s $output/$sample.variants.chr${chr}.raw.snv.vcf.multi.vcf ]]
				do
					echo " waiting for gatk and snvnix to complete to complete "
					sleep 2m	
				done
				### merge snvs and indels to give on vcf
				in="$output/$sample.variants.chr${chr}.raw.snv.vcf $output/$sample.variants.chr${chr}.raw.indel.vcf"
				$script_path/concatvcf.sh "$in" $output/$sample.variants.chr${chr}.raw.vcf $run_info yes
				in="$output/$sample.variants.chr${chr}.raw.snv.vcf.multi.vcf $output/$sample.variants.chr${chr}.raw.indel.vcf.multi.vcf"
				$script_path/concatvcf.sh "$in" $output/$sample.variants.chr${chr}.raw.multi.vcf $run_info yes
				$script_path/annotate_vcf.sh $output/$sample.variants.chr${chr}.raw.multi.vcf $run_info "$bam" 
				$script_path/annotate_vcf.sh $output/$sample.variants.chr${chr}.raw.vcf $run_info "$bam"
				
			fi
		elif [ $SNV_caller == "BEAUTY_EXOME" ]
		then
			## call indeles using GATK
			param="-L chr${chr}"
			bam="-I $output/$sample.chr${chr}-sorted.bam"
			$script_path/unifiedgenotyper.sh "$bam" $output/$sample.variants.chr${chr}.raw.gatk.vcf BOTH "$param" EMIT_VARIANTS_ONLY $run_info &
			### call snvs using snvmix
			$script_path/snvmix2.sh $sample "$bam" $output/$sample.variants.chr${chr}.raw.snvmix.vcf target "$param" $run_info &
			#UNION    
			while [[ ! -s $output/$sample.variants.chr${chr}.raw.gatk.vcf || ! -s $output/$sample.variants.chr${chr}.raw.snvmix.vcf  || ! -s $output/$sample.variants.chr${chr}.raw.gatk.vcf.multi.vcf || ! -s $output/$sample.variants.chr${chr}.raw.snvmix.vcf.multi.vcf ]]
			do
				echo " waiting for gatk and snvnix to complete to complete "
				sleep 2m	
			done
			input_var=""
			input_var="-V:GATK $output/$sample.variants.chr${chr}.raw.gatk.vcf -V:SNVMix $output/$sample.variants.chr${chr}.raw.snvmix.vcf -priority GATK,SNVMix"
			#Combine Variants
			$script_path/combinevcf.sh "$input_var" ${output}/$sample.variants.chr${chr}.raw.vcf $run_info yes
			$script_path/annotate_vcf.sh $output/$sample.variants.chr${chr}.raw.vcf $run_info "$bam" 
			input_var=""
			input_var="-V:GATK $output/$sample.variants.chr${chr}.raw.gatk.vcf.multi.vcf -V:SNVMix $output/$sample.variants.chr${chr}.raw.snvmix.vcf.multi.vcf -priority GATK,SNVMix"
			$script_path/combinevcf.sh "$input_var" ${output}/$sample.variants.chr${chr}.raw.multi.vcf $run_info yes
			$script_path/annotate_vcf.sh $output/$sample.variants.chr${chr}.raw.multi.vcf $run_info "$bam" 	
		fi
	else
		## assuming that normal is the first column/sample
		normal=${sampleArray[1]}.chr$chr-sorted.bam
		inputfiles="-I $output/$normal ";
		param="-L chr${chr}"
		if [ $somatic_caller == "BEAUTY_EXOME" ]
		then
			$script_path/snvmix2.sh ${sampleArray[1]} "$inputfiles" $output/${sampleArray[1]}.variants.chr${chr}.raw.snvmix.vcf target "$param" $run_info &
			in="-V $output/${sampleArray[1]}.variants.chr${chr}.raw.snvmix.vcf "
			in_mu="-V $output/${sampleArray[1]}.variants.chr${chr}.raw.snvmix.vcf.multi.vcf "
		fi
		for i in $(seq 2 ${#sampleArray[@]})
		do
			sample=${sampleArray[$i]}
			tumor=$sample.chr$chr-sorted.bam
			snv=$sample.chr$chr.snv.output
			inputfiles=$inputfiles" -I $output/$tumor"
			##run somatic caller
			if [ $somatic_caller == "JOINTSNVMIX" ]
			then
				$script_path/Jointsnvmix.sh $output/$normal $output/$tumor $output $chr $sample ${sampleArray[1]} $sample.chr$chr.snv.vcf $run_info
			elif [ $somatic_caller == "SOMATICSNIPER" ]
			then
				$script_path/somaticsnipper.sh $output/$normal $output/$tumor $output $chr $sample ${sampleArray[1]} $sample.chr$chr.snv.vcf $run_info
			elif [ $somatic_caller == "MUTECT" ]
			then
				$script_path/mutect.sh $output/$normal $output/$tumor $output $chr $sample ${sampleArray[1]} $sample.chr$chr.snv.vcf $run_info
			elif [ $somatic_caller == "BEAUTY_EXOME" ]
			then
				$script_path/mutect.sh $output/$normal $output/$tumor $output $chr $sample ${sampleArray[1]} $sample.chr$chr.snv.mutect.vcf $run_info
				$script_path/Jointsnvmix.sh $output/$normal $output/$tumor $output $chr $sample ${sampleArray[1]} $sample.chr$chr.snv.jsm.vcf $run_info & 
				$script_path/somaticsnipper.sh $output/$normal $output/$tumor $output $chr $sample ${sampleArray[1]} $sample.chr$chr.snv.ss.vcf $run_info & 
				while [[ ! -s $output/$sample.chr$chr.snv.jsm.vcf ||  ! -s $output/$sample.chr$chr.snv.ss.vcf ]]
				do
					echo " waiting for jointsnvmix somaticsniper to complete "
					sleep 2m	
				done
				#Combine vcf's into one VCF
				input_var="-V:MuTect $output/$sample.chr$chr.snv.mutect.vcf -V:JSM $output/$sample.chr$chr.snv.jsm.vcf -V:SomSniper $output/$sample.chr$chr.snv.ss.vcf -priority SomSniper,JSM,MuTect"
				#Combine Variants
				echo -e "starting to combine variants\n\n"
				$script_path/combinevcf.sh "$input_var" $output/$sample.chr${chr}.snv.vcf $run_info no
				#Add back the annotations that were lost (i.e. SomaticScore JSM_Prob
				### annotate vcfs
				in_bam="-I $input/$bam -resource:ss $output/$sample.chr$chr.snv.ss.vcf -resource:jsm $output/$sample.chr$chr.snv.jsm.vcf -resource:mutect $output/$sample.chr$chr.snv.mutect.vcf -E ss.SS -E ss.SSC -E jsm.PGERM -E jsm.PHETMUT -E jsm.PHOMMUT -E jsm.PLOH -E jsm.PPS -E jsm.PSOM -E mutect.MUTX -E mutect.POW -E mutect.MUTX_LOD"
				$script_path/annotate_vcf.sh $output/$sample.chr${chr}.snv.vcf $run_info "$in_bam"
				rm $output/$sample.chr$chr.snv.ss.vcf $output/$sample.chr$chr.snv.jsm.vcf $output/$sample.chr$chr.snv.mutect.vcf
				rm $output/$sample.chr$chr.snv.ss.vcf.idx $output/$sample.chr$chr.snv.jsm.vcf.idx $output/$sample.chr$chr.snv.mutect.vcf.idx
				input_var="-V:MuTect $output/$sample.chr$chr.snv.mutect.vcf.multi.vcf -V:JSM $output/$sample.chr$chr.snv.jsm.vcf.multi.vcf -V:SomSniper $output/$sample.chr$chr.snv.ss.vcf.multi.vcf -priority SomSniper,JSM,MuTect"
				#Combine Variants
				echo -e "starting to combine variants\n\n"
				$script_path/combinevcf.sh "$input_var" $output/$sample.chr${chr}.snv.vcf.multi.vcf $run_info no
				#Add back the annotations that were lost (i.e. SomaticScore JSM_Prob
				### annotate vcfs
				in_bam="-I $input/$bam -resource:ss $output/$sample.chr$chr.snv.ss.vcf.multi.vcf -resource:jsm $output/$sample.chr$chr.snv.jsm.vcf.multi.vcf -resource:mutect $output/$sample.chr$chr.snv.mutect.vcf.multi.vcf -E ss.SS -E ss.SSC -E jsm.PGERM -E jsm.PHETMUT -E jsm.PHOMMUT -E jsm.PLOH -E jsm.PPS -E jsm.PSOM -E mutect.MUTX -E mutect.POW -E mutect.MUTX_LOD"
				$script_path/annotate_vcf.sh $output/$sample.chr${chr}.snv.vcf.multi.vcf $run_info "$in_bam"
				rm $output/$sample.chr$chr.snv.ss.vcf.multi.vcf $output/$sample.chr$chr.snv.jsm.vcf.multi.vcf $output/$sample.chr$chr.snv.mutect.vcf.multi.vcf
				rm $output/$sample.chr$chr.snv.ss.vcf.multi.vcf.idx $output/$sample.chr$chr.snv.jsm.vcf.multi.vcf.idx $output/$sample.chr$chr.snv.mutect.vcf.multi.vcf.idx
				#remove old files
			fi
			### annotate vcfs if they haven't been already
			in_bam="-I $input/$bam"
			if [ $somatic_caller != "BEAUTY_EXOME" ]
			then
				$script_path/annotate_vcf.sh $output/$sample.chr$chr.snv.vcf $run_info "$in_bam"
			fi
			### somatic indel calling
			$script_path/somaticindel.sh $output/$tumor $output/$normal $chr $sample $output $sample.chr$chr.indel.vcf $run_info
			### annoatte vcfs
			$script_path/annotate_vcf.sh $output/$sample.chr$chr.indel.vcf $run_info "$in_bam"
		done
		## Germline calling
		param="-L chr${chr}"
		bam="-I $input/$bam"
		$script_path/unifiedgenotyper.sh "$bam" $output/variants.chr${chr}.raw.gatk.vcf BOTH "$param" EMIT_VARIANTS_ONLY $run_info
		### call snvs using snvmix
		if [ $somatic_caller == "BEAUTY_EXOME" ]
		then
			for i in $(seq 2 ${#sampleArray[@]})
			do
				sample=${sampleArray[$i]}
				tumor=$sample.chr$chr-sorted.bam
				in=$in"-V $output/$sample.variants.chr${chr}.raw.snvmix.vcf "
				in_mu=$in_mu"-V $output/$sample.variants.chr${chr}.raw.snvmix.vcf.multi.vcf "
				$script_path/snvmix2.sh $sample $output/$tumor $output/$sample.variants.chr${chr}.raw.snvmix.vcf target "$param" $run_info
			done
			while [[ ! -s $output/${sampleArray[1]}.variants.chr${chr}.raw.snvmix.vcf ]]
			do
				echo "waiting for vnsmix 2 complete for normal sample "
				sleep 2m	
			done
			$script_path/combinevcf.sh "$in" ${output}/variants.chr${chr}.raw.snvmix.vcf $run_info yes
			$script_path/combinevcf.sh "$in_mu" ${output}/variants.chr${chr}.raw.snvmix.vcf.multi.vcf $run_info yes
			#UNION
			input_var="-V:GATK $output/variants.chr${chr}.raw.gatk.vcf -V:SNVMIX ${output}/variants.chr${chr}.raw.snvmix.vcf -priority GATK,SNVMIX"
			#Combine Variants
			$script_path/combinevcf.sh "$input_var" ${output}/variants.chr${chr}.raw.vcf $run_info yes
			$script_path/annotate_vcf.sh $output/variants.chr${chr}.raw.vcf $run_info "$bam"
			input_var="-V:GATK $output/variants.chr${chr}.raw.gatk.vcf.multi.vcf -V:SNVMIX ${output}/variants.chr${chr}.raw.snvmix.vcf.multi.vcf -priority GATK,SNVMIX"
			#Combine Variants
			$script_path/combinevcf.sh "$input_var" ${output}/variants.chr${chr}.raw.multi.vcf $run_info yes
			$script_path/annotate_vcf.sh $output/variants.chr${chr}.raw.multi.vcf $run_info "$bam"
		else
			mv $output/variants.chr${chr}.raw.gatk.vcf ${output}/variants.chr${chr}.raw.vcf
			mv $output/variants.chr${chr}.raw.gatk.vcf.multi.vcf $output/variants.chr${chr}.raw.multi.vcf
			$script_path/annotate_vcf.sh ${output}/variants.chr${chr}.raw.vcf $run_info "$bam"
			$script_path/annotate_vcf.sh $output/variants.chr${chr}.raw.multi.vcf $run_info "$bam"
		fi		
	fi
	if [ $tool == "exome" ]
	then
		if [ -f $output/chr$chr.target.bed ]
		then
			rm $output/chr$chr.target.bed
		fi
	fi
	### after this we get multiple indels and snp files which need to be merged for Multi samples but just filter for
	if [ ${#sampleArray[@]} -gt 1 ]
	then
		### INDEL
		input_var=""
		multi_var=""
		for i in $(seq 2 ${#sampleArray[@]})
		do
			sample=${sampleArray[$i]}
			indel=$sample.chr$chr.indel.vcf
			multi=$sample.chr$chr.indel.vcf.multi.vcf
			input_var="${input_var} -V $output/$indel"
			multi_var="${multi_var} -V $output/$multi"
		done
		$script_path/combinevcf.sh "$input_var" $output/MergeAllSamples.chr$chr.Indels.raw.vcf $run_info yes
		$script_path/combinevcf.sh "$multi_var" $output/MergeAllSamples.chr$chr.Indels.raw.multi.vcf $run_info yes

		##Merge SNVs
		input_var=""
		multi_var=""
		for i in $(seq 2 ${#sampleArray[@]})
		do
			sample=${sampleArray[$i]}
			snv=$sample.chr$chr.snv.vcf
			multi=$sample.chr$chr.snv.vcf.multi.vcf
			input_var="${input_var} -V $output/$snv"
			multi_var="${multi_var} -V $output/$multi"
		done
		$script_path/combinevcf.sh "$input_var" $output/MergeAllSamples.chr$chr.snvs.raw.vcf $run_info yes
		$script_path/combinevcf.sh "$multi_var" $output/MergeAllSamples.chr$chr.snvs.raw.multi.vcf $run_info yes
		## combine both snv and indel
		in="$output/MergeAllSamples.chr$chr.snvs.raw.vcf $output/MergeAllSamples.chr$chr.Indels.raw.vcf"
		$script_path/concatvcf.sh "$in" $output/MergeAllSamples.chr$chr.raw.vcf $run_info yes
		in="$output/MergeAllSamples.chr$chr.snvs.raw.multi.vcf $output/MergeAllSamples.chr$chr.Indels.raw.multi.vcf"
		$script_path/concatvcf.sh "$in" $output/MergeAllSamples.chr$chr.raw.multi.vcf $run_info yes
	fi

	## remove files and add ED blat field
	if [ ${#sampleArray[@]} == 1 ]
	then
		### add the ED column
		if [[ ! -s $output/${sampleArray[1]}.variants.chr$chr.raw.vcf  || ! -s $output/${sampleArray[1]}.variants.chr$chr.raw.multi.vcf ]]
		then
			echo "ERROR : variant calling failed for ${sampleArray[1]} in variantss.h script"	
			exit 1;
		fi
		cat $output/${sampleArray[1]}.variants.chr$chr.raw.vcf | awk '$0 !~ /^#/ && $5 ~ /,/' > $output/${sampleArray[1]}.variants.chr$chr.raw.vcf.multi
        cat $output/${sampleArray[1]}.variants.chr$chr.raw.vcf | awk '$0 ~ /^#/ || $5 !~ /,/' > $output/${sampleArray[1]}.variants.chr$chr.raw.vcf.temp
        mv $output/${sampleArray[1]}.variants.chr$chr.raw.vcf.temp $output/${sampleArray[1]}.variants.chr$chr.raw.vcf
        cat $output/${sampleArray[1]}.variants.chr$chr.raw.multi.vcf $output/${sampleArray[1]}.variants.chr$chr.raw.vcf.multi > $output/${sampleArray[1]}.variants.chr$chr.raw.multi.vcf.temp
        $script_path/vcfsort.pl $ref.fai $output/${sampleArray[1]}.variants.chr$chr.raw.multi.vcf.temp > $output/${sampleArray[1]}.variants.chr$chr.raw.multi.vcf
        rm $output/${sampleArray[1]}.variants.chr$chr.raw.multi.vcf.temp $output/${sampleArray[1]}.variants.chr$chr.raw.vcf.multi
        $script_path/vcf_blat_verify.pl -i $output/${sampleArray[1]}.variants.chr$chr.raw.vcf -o $output/${sampleArray[1]}.variants.chr$chr.raw.vcf.temp -r $ref -b $blat -sam $samtools -br $blat_ref $blat_params
		mv $output/${sampleArray[1]}.variants.chr$chr.raw.vcf.temp $output/${sampleArray[1]}.variants.chr$chr.raw.vcf
		$script_path/vcf_blat_verify.pl -i $output/${sampleArray[1]}.variants.chr$chr.raw.multi.vcf -o $output/${sampleArray[1]}.variants.chr$chr.raw.multi.vcf.temp -r $ref -b $blat -sam $samtools -br $blat_ref $blat_params 
		mv $output/${sampleArray[1]}.variants.chr$chr.raw.multi.vcf.temp $output/${sampleArray[1]}.variants.chr$chr.raw.multi.vcf
		bam=chr${chr}.cleaned.bam
		rm $output/$bam.$chr.bam
		rm $output/$bam.$chr.bam.bai
		rm $output/${sampleArray[1]}.chr$chr.bam
		rm $output/${sampleArray[1]}.chr$chr.bam.bai
		rm $output/${sampleArray[1]}.chr$chr-sorted.bam
		rm $output/${sampleArray[1]}.chr$chr-sorted.bam.bai
		### update the file size
		$script_path/filesize.sh VariantCalling ${sampleArray[1]} $output ${sampleArray[1]}.variants.chr$chr.raw.vcf $JOB_ID $run_info
	fi
	if [ ${#sampleArray[@]} -gt 1 ]
	then
		if [[ ! -s $output/variants.chr$chr.raw.vcf || ! -s $output/variants.chr$chr.raw.multi.vcf || ! -s $output/MergeAllSamples.chr$chr.raw.vcf || ! -s $output/MergeAllSamples.chr$chr.raw.multi.vcf ]]
		then
			echo "ERROR : variant calling failed for pair in variants.sh script"	
			exit 1;
		fi
		cat $output/variants.chr$chr.raw.vcf |  awk '$0 !~ /^#/ && $5 ~ /,/' > $output/variants.chr$chr.raw.vcf.multi
		cat $output/variants.chr$chr.raw.vcf |  awk '$0 ~ /^#/ || $5 !~ /,/' > $output/variants.chr$chr.raw.vcf.temp
		mv $output/variants.chr$chr.raw.vcf.temp $output/variants.chr$chr.raw.vcf
		cat $output/variants.chr$chr.raw.multi.vcf $output/variants.chr$chr.raw.vcf.multi > $output/variants.chr$chr.raw.multi.vcf.temp
		$script_path/vcfsort.pl $ref.fai $output/variants.chr$chr.raw.multi.vcf.temp > $output/variants.chr$chr.raw.multi.vcf
		rm $output/variants.chr$chr.raw.multi.vcf.temp $output/variants.chr$chr.raw.vcf.multi
		$script_path/vcf_blat_verify.pl -i $output/variants.chr$chr.raw.vcf -o $output/variants.chr$chr.raw.vcf.temp -r $ref -b $blat -sam $samtools -br $blat_ref $blat_params
		mv $output/variants.chr$chr.raw.vcf.temp $output/variants.chr$chr.raw.vcf
		$script_path/vcf_blat_verify.pl -i $output/variants.chr$chr.raw.multi.vcf -o $output/variants.chr$chr.raw.multi.vcf.temp -r $ref -b $blat -sam $samtools -br $blat_ref $blat_params
		mv $output/variants.chr$chr.raw.multi.vcf.temp $output/variants.chr$chr.raw.multi.vcf
		
		cat $output/MergeAllSamples.chr$chr.raw.vcf |  awk '$0 !~ /^#/ && $5 ~ /,/' > $output/MergeAllSamples.chr$chr.raw.vcf.multi
		cat $output/MergeAllSamples.chr$chr.raw.vcf |  awk '$0 ~ /^#/ || $5 !~ /,/' > $output/MergeAllSamples.chr$chr.raw.vcf.temp
		mv $output/MergeAllSamples.chr$chr.raw.vcf.temp $output/MergeAllSamples.chr$chr.raw.vcf
		cat $output/MergeAllSamples.chr$chr.raw.multi.vcf $output/MergeAllSamples.chr$chr.raw.vcf.multi > $output/MergeAllSamples.chr$chr.raw.multi.vcf.temp
		$script_path/vcfsort.pl $ref.fai $output/MergeAllSamples.chr$chr.raw.multi.vcf.temp > $output/MergeAllSamples.chr$chr.raw.multi.vcf
		rm $output/MergeAllSamples.chr$chr.raw.multi.vcf.temp $output/MergeAllSamples.chr$chr.raw.vcf.multi
		$script_path/vcf_blat_verify.pl -i $output/MergeAllSamples.chr$chr.raw.vcf -o $output/MergeAllSamples.chr$chr.raw.vcf.temp -r $ref -b $blat -sam $samtools -br $blat_ref $blat_params
		mv $output/MergeAllSamples.chr$chr.raw.vcf.temp $output/MergeAllSamples.chr$chr.raw.vcf
		$script_path/vcf_blat_verify.pl -i $output/MergeAllSamples.chr$chr.raw.multi.vcf -o $output/MergeAllSamples.chr$chr.raw.multi.vcf.temp -r $ref -b $blat -sam $samtools -br $blat_ref $blat_params
		mv $output/MergeAllSamples.chr$chr.raw.multi.vcf.temp $output/MergeAllSamples.chr$chr.raw.multi.vcf 
		for i in $(seq 1 ${#sampleArray[@]})
		do	
			rm $output/${sampleArray[$i]}.chr$chr.rg.bam
			rm $output/${sampleArray[$i]}.chr$chr.rg.bam.$chr.bam
			rm $output/${sampleArray[$i]}.chr$chr.rg.bam.$chr.bam.bai
			rm $output/${sampleArray[$i]}.chr$chr.bam
			rm $output/${sampleArray[$i]}.chr$chr.bam.bai
			rm $output/${sampleArray[$i]}.chr$chr-sorted.bam
			rm $output/${sampleArray[$i]}.chr$chr-sorted.bam.bai
		done
		$script_path/filesize.sh VariantCalling multi_sample $output variants.chr$chr.raw.vcf $JOB_ID $run_info
		$script_path/filesize.sh VariantCalling multi_sample $output MergeAllSamples.chr$chr.raw.vcf $JOB_ID $run_info
	fi
	## update dash board
	if [ $SGE_TASK_ID == 1 ]
	then
		for i in `echo $samples | tr ":" " "`
		do
			$script_path/dashboard.sh $i $run_info VariantCalling complete
		done
	fi
	echo `date`
fi
