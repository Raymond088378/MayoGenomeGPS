#!/bin/bash

if [ $# -le 4 ]
then
	echo -e "script to run variant calling on set of BAM files\nUsage: ./variants.sh \
</path/to/input dir> <samples ':' sep[normal:tumor1:tumor2:tumorN]> </path/to/output directory> \
<1 or 0 chopped or not > <path/to/run info> <SGE_TASK_ID(optional)>\
\nNOTE: first sample is considered as normal and others are considered as tumor[Assumption]";
	exit 1;
fi

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
somatic_calling=$( cat $tool_info | grep -w '^SOMATIC_CALLING' | cut -d '=' -f2 | tr "[a-z]" "[A-Z]" )
export PERL5LIB=$perllib:$PERL5LIB
export PATH=$java:$tabix/:$PATH

#### check and validate the bam file and let user to proceed after validation
bam=chr${chr}.cleaned.bam
$samtools/samtools view -H $input/$bam 1>$input/$bam.$chr.header 2> $input/$bam.$chr.fix.log
if [[ `cat $input/$bam.$chr.fix.log | wc -l` -gt 0 || `cat $input/$bam.$chr.header | wc -l` -le 0 ]]
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

### update dashborad
if [ $SGE_TASK_ID == 1 ]
then
	for i in `echo $samples | tr ":" " "`
		do
			$script_path/dashboard.sh $i $run_info VariantCalling started
		done
fi

### fill array with sample names
i=1
for sample in $sampleNames
do
	sampleArray[$i]=$sample
	let i=i+1
done

### DEAD CODE
### num_samples=`echo $samples | tr ":" " "`

mkdir -p $output/temp

### if multiple samples, validate
### if multiple and somatic calling, then split BAMs using read group and validate 
if [ ${#sampleArray[@]} -gt 1 ]
then
	#bams are splitted using read group information
	if [ $somatic_calling == "YES" ]
	then
		for i in $(seq 1 ${#sampleArray[@]})
		do
			sample=${sampleArray[$i]}
			#removing the header for the extra samples from the BAM
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
	fi
else
	sample=${sampleArray[1]}
		$script_path/samplecheckBAM.sh $input $bam $output $run_info $sample $chopped $chr
	fi

	inputfiles=""
if [ $tool == "exome" ]
then
	cat $TargetKit | grep -w chr$chr > $output/chr$chr.target.bed
fi

### OUTPUT (SINGLE SAMPLE:GATK:EXOME:ALLSITES) : $output/$sample.variants.chr${chr}.raw.all.vcf.bgz  (no ./. calls)
### 											 $output/$sample.variants.chr${chr}.raw.vcf
### OUTPUT (SINGLE SAMPLE:GATK:EXOME:KNOWN SITES) : $output/$sample.variants.chr${chr}.raw.vcf
### OUTPUT (SINGLE SAMPLE:GATK:WHOLEGENOME:KNOWN SITES) : $output/$sample.variants.chr${chr}.raw.vcf 
			


### If we only have one sample
if [ ${#sampleArray[@]} == 1 ] 
then
	sample=${sampleArray[1]}
	if [ $SNV_caller == "GATK" ]
	then
		## call variants using gatk UnifiedGenotyper module
		if [[ $all_sites == "YES"  && $tool == "exome" ]]
		then
			len=`cat $output/chr$chr.target.bed |wc -l`
			if [ $len -gt 0 ]
			then    
				param="-L $output/chr$chr.target.bed"
			else
				param="-L chr$chr"
			fi
			bam="-I $output/$sample.chr${chr}-sorted.bam"
			$script_path/unifiedgenotyper.sh "$bam" $output/$sample.variants.chr${chr}.raw.all.vcf BOTH "$param" EMIT_ALL_SITES $run_info
			
			### add phase by transmission if pedigree information provided.
			if [ $ped != "NA" ]
			then
				$script_path/phaseByTransmission.sh $output/$sample.variants.chr${chr}.raw.all.vcf $output/$sample.variants.chr${chr}.raw.all.pbt.vcf $run_info
			fi
			
			### filter raw.all.vcf file to keep only the variant calls
			cat $output/$sample.variants.chr${chr}.raw.all.vcf | awk '$5 != "." || $0 ~ /^#/' | grep -v "\./\."  | grep -v "0\/0" > $output/$sample.variants.chr${chr}.raw.vcf
			sed '/^$/d' $output/$sample.variants.chr${chr}.raw.vcf > $output/$sample.variants.chr${chr}.raw.vcf.temp
			mv $output/$sample.variants.chr${chr}.raw.vcf.temp $output/$sample.variants.chr${chr}.raw.vcf
			
			### prepare the file for backfilling, delete all lines containing ./. 
			cat $output/$sample.variants.chr${chr}.raw.all.vcf | grep -v "\./\." > $output/$sample.variants.chr${chr}.raw.all.vcf.temp
			mv $output/$sample.variants.chr${chr}.raw.all.vcf.temp $output/$sample.variants.chr${chr}.raw.all.vcf
			
			### remove all empty lines
			sed '/^$/d' $output/$sample.variants.chr${chr}.raw.all.vcf > $output/$sample.variants.chr${chr}.raw.all.vcf.temp
			mv $output/$sample.variants.chr${chr}.raw.all.vcf.temp $output/$sample.variants.chr${chr}.raw.all.vcf
			
			### call tabix and zip up the vcf
			$tabix/bgzip $output/$sample.variants.chr${chr}.raw.all.vcf
			
		elif [ $tool == "exome" ] ### all_sites == "NO"
		then
			len=`cat $output/chr$chr.target.bed |wc -l`
			if [ $len -gt 0 ]
			then    
				param="-L $output/chr$chr.target.bed"
			else
				param="-L chr$chr"
			fi
			bam="-I $output/$sample.chr${chr}-sorted.bam"
			$script_path/unifiedgenotyper.sh "$bam" $output/$sample.variants.chr${chr}.raw.vcf BOTH "$param" EMIT_VARIANTS_ONLY $run_info
			
			
			## add phase by transmission if pedigree information provided.
			if [ $ped != "NA" ]
			then
				### FIXED: THIS REFERED TO A FILE THAT SHOULD NOT EXIST AT THIS POINT (.all.vcf) CR, 2/14/2013
				$script_path/phaseByTransmission.sh $output/$sample.variants.chr${chr}.raw.vcf $output/$sample.variants.chr${chr}.raw.pbt.vcf $run_info
			fi
			$script_path/annotate_vcf.sh $output/$sample.variants.chr${chr}.raw.vcf $run_info "$bam"
		else 
			### tool == whole_genome ASSUMED all_sites=="NO"
			param="-L chr$chr"
			bam="-I $output/$sample.chr${chr}-sorted.bam"
			$script_path/unifiedgenotyper.sh "$bam" $output/$sample.variants.chr${chr}.raw.vcf BOTH "$param" EMIT_VARIANTS_ONLY $run_info
			## add phase by transmission if pedigree information provided.
			if [ $ped != "NA" ]
			then
				### FIXED: THIS REFERED TO A FILE THAT SHOULD NOT EXIST AT THIS POINT (.all.vcf) CR, 2/14/2013
				$script_path/phaseByTransmission.sh $output/$sample.variants.chr${chr}.raw.vcf $output/$sample.variants.chr${chr}.raw.pbt.vcf $run_info
			fi
			
			$script_path/annotate_vcf.sh $output/$sample.variants.chr${chr}.raw.vcf $run_info "$bam"
		fi
	elif [ $SNV_caller == "SNVMIX" ]
	then
		### call indels using GATK
		if [[ $all_sites == "YES"  && $tool == "exome" ]]
		then
			len=`cat $output/chr$chr.target.bed |wc -l`
			if [ $len -gt 0 ]
			then    
				param="-L $output/chr$chr.target.bed"
			else
				param="-L chr$chr"
			fi
			## call indels
			bam="-I $output/$sample.chr${chr}-sorted.bam"
			$script_path/unifiedgenotyper.sh "$bam" $output/$sample.variants.chr${chr}.raw.indel.all.vcf INDEL "$param" EMIT_ALL_SITES $run_info
			### call snvs using SNVmix
			$script_path/snvmix2.sh $sample "$bam" $output/$sample.variants.chr${chr}.raw.snv.all.vcf all "$param" $run_info
			### annoatte vcf
			$script_path/annotate_vcf.sh $output/$sample.variants.chr${chr}.raw.snv.all.vcf $run_info "$bam"
			### merge snvs and indels to give .raw.all.vcf
			in="$output/$sample.variants.chr${chr}.raw.snv.all.vcf $output/$sample.variants.chr${chr}.raw.indel.all.vcf "
			$script_path/concatvcf.sh "$in" $output/$sample.variants.chr${chr}.raw.all.vcf $run_info yes
			
			### Remove uncalled variants to give .raw.vcf
			cat $output/$sample.variants.chr${chr}.raw.all.vcf | awk '$5 != "." || $0 ~ /^#/' | grep -v "\./\."  | grep -v "0\/0" > $output/$sample.variants.chr${chr}.raw.vcf
			$tabix/bgzip $output/$sample.variants.chr${chr}.raw.all.vcf
			
		elif [ $tool == "exome" ]
		then
			len=`cat $output/chr$chr.target.bed |wc -l`
			if [ $len -gt 0 ]
			then    
				param="-L $output/chr$chr.target.bed"
			else
				param="-L chr$chr"
			fi
			bam="-I $output/$sample.chr${chr}-sorted.bam"
			$script_path/unifiedgenotyper.sh "$bam" $output/$sample.variants.chr${chr}.raw.indel.vcf INDEL "$param" EMIT_VARIANTS_ONLY $run_info &
			### call snvs using snvmix
			$script_path/snvmix2.sh $sample "$bam" $output/$sample.variants.chr${chr}.raw.snv.vcf target "$param" $run_info &
			while [[ ! -s $output/$sample.variants.chr${chr}.raw.indel.vcf || ! -s $output/$sample.variants.chr${chr}.raw.snv.vcf ]]
			do
				echo " waiting for gatk and snvnix to complete to complete "
				sleep 2m	
			done
			### merge snvs and indels to give on vcf
			in="$output/$sample.variants.chr${chr}.raw.snv.vcf $output/$sample.variants.chr${chr}.raw.indel.vcf"
			$script_path/concatvcf.sh "$in" $output/$sample.variants.chr${chr}.raw.vcf $run_info yes
			$script_path/annotate_vcf.sh $output/$sample.variants.chr${chr}.raw.vcf $run_info "$bam"
		else
			param="-L chr$chr"
			bam="-I $output/$sample.chr${chr}-sorted.bam"
			$script_path/unifiedgenotyper.sh "$bam" $output/$sample.variants.chr${chr}.raw.indel.vcf INDEL "$param" EMIT_VARIANTS_ONLY $run_info &
			### call snvs using snvmix
			$script_path/snvmix2.sh $sample "$bam" $output/$sample.variants.chr${chr}.raw.snv.vcf target "$param" $run_info &
			while [[ ! -s $output/$sample.variants.chr${chr}.raw.indel.vcf  || ! -s $output/$sample.variants.chr${chr}.raw.snv.vcf ]]
			do
				echo " waiting for gatk and snvnix to complete to complete "
				sleep 2m	
			done
			### merge snvs and indels to give on vcf
			in="$output/$sample.variants.chr${chr}.raw.snv.vcf $output/$sample.variants.chr${chr}.raw.indel.vcf"
			$script_path/concatvcf.sh "$in" $output/$sample.variants.chr${chr}.raw.vcf $run_info yes
			$script_path/annotate_vcf.sh $output/$sample.variants.chr${chr}.raw.vcf $run_info "$bam"
		fi
		### SNVMix on exome or whole genome
	elif [ $SNV_caller == "BEAUTY_EXOME" ]
	then
		if [ $tool == "exome" ]
		then    
			len=`cat $output/chr$chr.target.bed |wc -l`
			if [ $len -gt 0 ]
			then    
				param="-L $output/chr$chr.target.bed"
			else
				param="-L chr$chr"
			fi
			bam="-I $output/$sample.chr${chr}-sorted.bam"
			$script_path/unifiedgenotyper.sh "$bam" $output/$sample.variants.chr${chr}.raw.gatk.vcf BOTH "$param" EMIT_VARIANTS_ONLY $run_info &
			### call snvs using snvmix
			$script_path/snvmix2.sh $sample "$bam" $output/$sample.variants.chr${chr}.raw.snvmix.vcf target "$param" $run_info &
			
			### TODO: Change to pid/wait instead of unbound loop
			while [[ ! -s $output/$sample.variants.chr${chr}.raw.gatk.vcf || ! -s $output/$sample.variants.chr${chr}.raw.snvmix.vcf  ]]
			do
				echo " waiting for gatk and snvnix to complete to complete "
				sleep 2m	
			done
			
			input_var=""
			input_var="-V:GATK $output/$sample.variants.chr${chr}.raw.gatk.vcf -V:SNVMix $output/$sample.variants.chr${chr}.raw.snvmix.vcf -priority GATK,SNVMix"
			#Combine Variants
			$script_path/combinevcf.sh "$input_var" ${output}/$sample.variants.chr${chr}.raw.vcf $run_info yes
			$script_path/annotate_vcf.sh $output/$sample.variants.chr${chr}.raw.vcf $run_info "$bam" 
		else
			param="-L chr$chr"
			bam="-I $output/$sample.chr${chr}-sorted.bam"
			$script_path/unifiedgenotyper.sh "$bam" $output/$sample.variants.chr${chr}.raw.gatk.vcf BOTH "$param" EMIT_VARIANTS_ONLY $run_info &
			### call snvs using snvmix
			$script_path/snvmix2.sh $sample "$bam" $output/$sample.variants.chr${chr}.raw.snvmix.vcf target "$param" $run_info &
			
			### TODO: Change to pid / wait instead of unbound loop
			while [[ ! -s $output/$sample.variants.chr${chr}.raw.gatk.vcf || ! -s $output/$sample.variants.chr${chr}.raw.snvmix.vcf ]]
			do
				echo " waiting for gatk and snvnix to complete to complete "
				sleep 2m	
			done
			
			#UNION    
			input_var=""
			input_var="-V:GATK $output/$sample.variants.chr${chr}.raw.gatk.vcf -V:SNVMix $output/$sample.variants.chr${chr}.raw.snvmix.vcf -priority GATK,SNVMix"
			#Combine Variants
			$script_path/combinevcf.sh "$input_var" ${output}/$sample.variants.chr${chr}.raw.vcf $run_info yes
			$script_path/annotate_vcf.sh $output/$sample.variants.chr${chr}.raw.vcf $run_info "$bam"   
		fi ### BEAUTY exome or whole genome
	fi ### SNV Callers 
else
	### Multiple samples, assumption that normal is the first column/sample
	normal=${sampleArray[1]}.chr$chr-sorted.bam
	inputfiles="-I $output/$normal ";
	
	if [ $tool == "exome" ]
	then
		len=`cat $output/chr$chr.target.bed |wc -l`
		if [ $len -gt 0 ]
		then    
			param="-L $output/chr$chr.target.bed"
		else
			param="-L chr$chr"
		fi
	else                
		param="-L chr${chr}"
	fi
	
	if [ $somatic_caller == "BEAUTY_EXOME" ]
	then
		### Start snvmix calling on normal sample 0xDEADBEEF (get pid here, wait later)
		$script_path/snvmix2.sh ${sampleArray[1]} "$inputfiles" $output/${sampleArray[1]}.variants.chr${chr}.raw.snvmix.vcf target "$param" $run_info &
		in="-V $output/${sampleArray[1]}.variants.chr${chr}.raw.snvmix.vcf "
	fi
	
	if [ $somatic_calling == "YES" ]
	then
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
				
				### Start mutect, jointsnv, somatic sniper calling on samples 2..N
				$script_path/mutect.sh $output/$normal $output/$tumor $output $chr $sample ${sampleArray[1]} $sample.chr$chr.snv.mutect.vcf $run_info
				$script_path/Jointsnvmix.sh $output/$normal $output/$tumor $output $chr $sample ${sampleArray[1]} $sample.chr$chr.snv.jsm.vcf $run_info & 
				$script_path/somaticsnipper.sh $output/$normal $output/$tumor $output $chr $sample ${sampleArray[1]} $sample.chr$chr.snv.ss.vcf $run_info &
				
				
				### TODO: Change this to a PID / wait statement instead of polling for file creation
				while [[ ! -s $output/$sample.chr$chr.snv.jsm.vcf ||  ! -s $output/$sample.chr$chr.snv.ss.vcf ]]
				do
					echo " waiting for jointsnvmix somaticsniper to complete "
					sleep 2m	
				done
				
				### combine INPUT .snv.mutect.vcf, .snv.jsm.vcf, .snv.ss.vcf
				input_var="-V:MuTect $output/$sample.chr$chr.snv.mutect.vcf -V:JSM $output/$sample.chr$chr.snv.jsm.vcf -V:SomSniper $output/$sample.chr$chr.snv.ss.vcf -priority SomSniper,JSM,MuTect"
				echo -e "starting to combine variants\n\n"
				$script_path/combinevcf.sh "$input_var" $output/$sample.chr${chr}.snv.vcf $run_info yes
				### OUTPUT -- $output/$sample.chr${chr}.snv.vcf
			fi
			
			### annotate snvs (Beauty handles later)
			in_bam="-I $input/$bam"
			if [ $somatic_caller != "BEAUTY_EXOME" ]
			then
				$script_path/annotate_vcf.sh $output/$sample.chr$chr.snv.vcf $run_info "$in_bam"
			fi
			
			### somatic indel calling
			$script_path/somaticindel.sh $output/$tumor $output/$normal $chr "$param" $sample $output $sample.chr$chr.indel.vcf $run_info
			
			### annotate indels 
			$script_path/annotate_vcf.sh $output/$sample.chr$chr.indel.vcf $run_info "$in_bam"
			
			### OUTPUT -- $output/$sample.chr$chr.indel.vcf
		done
		### finished samples 2..N 
	fi
	### endif somatic calling

	##########################################
	### Germline Calling on Main .bam File ###
	### Multiple Samples Case              ###
	##########################################

	bam="-I $input/$bam"
	$script_path/unifiedgenotyper.sh "$bam" $output/variants.chr${chr}.raw.gatk.vcf BOTH "$param" EMIT_VARIANTS_ONLY $run_info
	
	#### $output/variants.chr${chr}.raw.gatk.vcf is later merged into ${output}/variants.chr${chr}.raw.vcf
	
	if [ $somatic_caller == "BEAUTY_EXOME" ]
	then
		### call snvs using snvmix on samples 2..N
		### could this be relocated into above loop?? 
		for i in $(seq 2 ${#sampleArray[@]})
		do
			sample=${sampleArray[$i]}
			tumor=$sample.chr$chr-sorted.bam
			in=$in"-V $output/$sample.variants.chr${chr}.raw.snvmix.vcf "
			$script_path/snvmix2.sh $sample $output/$tumor $output/$sample.variants.chr${chr}.raw.snvmix.vcf target "$param" $run_info
		done
		
		### TODO: Replace with PID/wait instead of unbounded loop see:0xDEADBEEF
		while [[ ! -s $output/${sampleArray[1]}.variants.chr${chr}.raw.snvmix.vcf ]]
		do
			echo "waiting for snvmix2 complete for normal sample "
			sleep 2m	
		done
		
		$script_path/combinevcf.sh "$in" ${output}/variants.chr${chr}.raw.snvmix.vcf $run_info yes
		
		#UNION - combine .gatk.vcf, .snvmix.vcf 
		input_var="-V:GATK $output/variants.chr${chr}.raw.gatk.vcf -V:SNVMIX ${output}/variants.chr${chr}.raw.snvmix.vcf -priority GATK,SNVMIX"
		#Combine Variants
		$script_path/combinevcf.sh "$input_var" ${output}/variants.chr${chr}.raw.vcf $run_info yes
		$script_path/annotate_vcf.sh $output/variants.chr${chr}.raw.vcf $run_info "$bam"
	else
		mv $output/variants.chr${chr}.raw.gatk.vcf ${output}/variants.chr${chr}.raw.vcf
		$script_path/annotate_vcf.sh ${output}/variants.chr${chr}.raw.vcf $run_info "$bam"
	fi		
	
fi ### Done with multiple sample part of script

### XXX here we should have $output/variants.chr${chr}.raw.vcf no matter what XXX

if [ $tool == "exome" ]
then
	if [ -f $output/chr$chr.target.bed ]
	then
		rm $output/chr$chr.target.bed
	fi
fi

### after this we get multiple indels and snp files which need to be merged and backfilled for somatic samples 
if [ ${#sampleArray[@]} -gt 1 ]
then
	if [ $somatic_calling == "YES" ]
	then
		input_var=""
		### Merge INDEL Files Across Samples 2..N
		for i in $(seq 2 ${#sampleArray[@]})
		do
			sample=${sampleArray[$i]}
			indel=$sample.chr$chr.indel.vcf
			input_var="${input_var} -V $output/$indel"
		done
		$script_path/combinevcf.sh "$input_var" $output/MergeAllSamples.chr$chr.Indels.raw.vcf $run_info yes

		###Merge SNVs Across Samples 2..N
		input_var=""
		for i in $(seq 2 ${#sampleArray[@]})
		do
			sample=${sampleArray[$i]}
			snv=$sample.chr$chr.snv.vcf
			input_var="${input_var} -V $output/$snv"
		done
		
		$script_path/combinevcf.sh "$input_var" $output/MergeAllSamples.chr$chr.snvs.raw.vcf $run_info yes
		
		### Perform Backfilling on SNVs in Samples 2..N
		$script_path/unifiedgenotyper_backfill.sh "-I $input/chr${chr}.cleaned.bam" $output/MergeAllSamples.chr$chr.snvs.raw.vcf BOTH EMIT_ALL_SITES $run_info
		
		## combine both snv and indel
		in="$output/MergeAllSamples.chr$chr.snvs.raw.vcf $output/MergeAllSamples.chr$chr.Indels.raw.vcf"
		$script_path/concatvcf.sh "$in" $output/MergeAllSamples.chr$chr.raw.vcf $run_info yes
	fi
fi

### Cleanup single sample
if [ ${#sampleArray[@]} == 1 ]
then
	### remove files and add ED blat field column
	if [[ ! -s $output/${sampleArray[1]}.variants.chr$chr.raw.vcf ]]
	then
		echo "ERROR : variant calling failed for ${sampleArray[1]} in variants.sh script"	
		exit 1;
	fi
	$script_path/vcf_blat_verify.pl -i $output/${sampleArray[1]}.variants.chr$chr.raw.vcf -o $output/${sampleArray[1]}.variants.chr$chr.raw.vcf.temp -r $ref -b $blat -sam $samtools -br $blat_ref $blat_params
	mv $output/${sampleArray[1]}.variants.chr$chr.raw.vcf.temp $output/${sampleArray[1]}.variants.chr$chr.raw.vcf
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

### Cleanup multiple samples
if [ ${#sampleArray[@]} -gt 1 ]
then
	if [ $somatic_calling == "NO" ]
	then
		if [ ! -s $output/variants.chr$chr.raw.vcf ]
		then
			echo "ERROR : variant calling failed for pair in variants.sh script"
		fi	
	else
		if [[ ! -s $output/variants.chr$chr.raw.vcf || ! -s $output/MergeAllSamples.chr$chr.raw.vcf ]]
		then
			echo "ERROR : variant calling failed for pair in variants.sh script"	
			exit 1;
		fi	
	fi
	$script_path/vcf_blat_verify.pl -i $output/variants.chr$chr.raw.vcf -o $output/variants.chr$chr.raw.vcf.temp -r $ref -b $blat -sam $samtools -br $blat_ref $blat_params
	mv $output/variants.chr$chr.raw.vcf.temp $output/variants.chr$chr.raw.vcf
	if [ $somatic_calling == "YES" ]
	then
		$script_path/vcf_blat_verify.pl -i $output/MergeAllSamples.chr$chr.raw.vcf -o $output/MergeAllSamples.chr$chr.raw.vcf.temp -r $ref -b $blat -sam $samtools -br $blat_ref $blat_params
		mv $output/MergeAllSamples.chr$chr.raw.vcf.temp $output/MergeAllSamples.chr$chr.raw.vcf
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
		$script_path/filesize.sh VariantCalling multi_sample $output MergeAllSamples.chr$chr.raw.vcf $JOB_ID $run_info
	fi
	$script_path/filesize.sh VariantCalling multi_sample $output variants.chr$chr.raw.vcf $JOB_ID $run_info
fi

### update dash board
if [ $SGE_TASK_ID == 1 ]
then
	for i in `echo $samples | tr ":" " "`
	do
		$script_path/dashboard.sh $i $run_info VariantCalling complete
	done
fi
echo `date`
