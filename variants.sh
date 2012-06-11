

if [ $# -le 4 ]
then
    echo -e "Usage:\n<input dir><samples ':' sep[normal:tumor1:tumor2:tumorN]> <output directory> <1 or 0 chopped or not ><run info>\nNOTE: first sample is considered as normal and others are considered as tumor[Assumption]\n";
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
    script_path=$( cat $tool_info | grep -w '^WHOLEGENOME_PATH' | cut -d '=' -f2 )
    analysis=$( cat $run_info | grep -w '^ANALYSIS' | cut -d '=' -f2 |tr "[A-Z]" "[a-z]")
    sampleNames=$( echo $samples | tr ":" "\n" )
    chr=$( cat $run_info | grep -w '^CHRINDEX' | cut -d '=' -f2 | tr ":" "\n" | head -n $SGE_TASK_ID | tail -n 1 )
    all_sites=$( cat $tool_info | grep -w '^EMIT_ALL_SITES' | cut -d '=' -f2 | tr "[a-z]" "[A-Z]" )
    depth_filter=$( cat $tool_info | grep -w '^DEPTH_FILTER' | cut -d '=' -f2)
    prob_filter=$( cat $tool_info | grep -w '^PROB_FILTER' | cut -d '=' -f2)
    only_ontarget=$( cat $tool_info | grep -w '^TARGETTED' | cut -d '=' -f2 | tr "[a-z]" "[A-Z]" )
    TargetKit=$( cat $tool_info | grep -w '^ONTARGET' | cut -d '=' -f2 )
    SNV_caller=$( cat $run_info | grep -w '^SNV_CALLER' | cut -d '=' -f2)
    somatic_caller=$( cat $run_info | grep -w '^SOMATIC_CALLER' | cut -d '=' -f2)
	ped=$( cat $tool_info | grep -w '^PEDIGREE' | cut -d '=' -f2 )
	javahome=$( cat $tool_info | grep -w '^JAVA_HOME' | cut -d '=' -f2 )
	
	export PERL5LIB=$PERL5LIB:$perllib
    export PATH=$tabix/:$PATH

	export JAVA_HOME=$javahome
	export PATH=$JAVA_HOME/bin:$PATH
	
	
    bam=chr${chr}.cleaned.bam

    ## update dashborad
    if [ $SGE_TASK_ID == 1 ]
    then
		$script_path/dashboard.sh $samples $run_info VariantCalling started
    fi

    i=1
    for sample in $sampleNames
    do
       sampleArray[$i]=$sample
       let i=i+1
    done
    num_samples=`echo $samples | tr ":" " "`

    ## if multiple samples then split using read group and Validate BAM other wise just check and validate BAM
    if [ ${#sampleArray[@]} -gt 1 ]
    then
        ## bams are splitted using read group information
        for i in $(seq 1 ${#sampleArray[@]})
        do
            sample=${sampleArray[$i]}
            ### removing the header for the extra samples from the BAM
            sam=`echo $samples | tr ":" "\n"| grep -v "$sample" | tr "\n" " "`
            gr=""
            for s in $sam
            do
                a="ID:$s|";
                gr="$gr $a"
            done
            gr=`echo $gr |  sed "s/|$//"`
            $samtools/samtools view -b -r $sample $input/$bam > $output/$sample.chr$chr.rg.bam
            $samtools/samtools view -H $output/$sample.chr$chr.rg.bam | grep -E -v "$gr" | $samtools/samtools reheader - $output/$sample.chr$chr.rg.bam > $output/$sample.chr$chr.rg.re.bam
            mv $output/$sample.chr$chr.rg.re.bam $output/$sample.chr$chr.rg.bam

            if [ ! -s $output/$sample.chr$chr.rg.bam ]
            then
                $script_path/errorlog.sh $output/$sample.chr$chr.rg.bam variants.sh ERROR "failed to create"
                exit 1
            fi
            $script_path/samplecheckBAM.sh $output $sample.chr$chr.rg.bam $output $run_info $sample $chopped $chr
        done
    else
        sample=${sampleArray[1]}
        if [ $SNV_caller == SNVMIX ]
        then
            $samtools/samtools pileup -s -f $ref $input/$bam > $input/chr${chr}.pileup
        fi
        $script_path/samplecheckBAM.sh $input $bam $output $run_info $sample $chopped $chr
    fi

    inputfiles=""
    if [ ${#sampleArray[@]} == 1 ]
    then
        sample=${sampleArray[1]}
        ### use GATK or SNVmix
        if [ $SNV_caller == "GATK" ]
        then
            ## call variants using gatk UnifiedGenotyper module
            if [[ $all_sites == "YES"  && $tool == "exome" ]]
            then
                cat $TargetKit | grep -w chr$chr > $output/$sample.$chr.target.bed
                len=`cat $output/$sample.$chr.target.bed |wc -l`
                if [[ $only_ontarget == "YES" && $len -gt 0 ]]
                then
                    param="-L $output/$sample.$chr.target.bed"
                else
                    param="-L chr${chr}"
                fi
                bam="-I $output/$sample.chr${chr}-sorted.bam"
                $script_path/unifiedgenotyper.sh "$bam" $output/$sample.variants.chr${chr}.raw.all.vcf BOTH "$param" EMIT_ALL_SITES $run_info
                rm $output/$sample.$chr.target.bed

				## add phase by transmission if pedigree information provided.
				if [$ped != "NA"]
				then
					$script_path/phaseByTransmission.sh $output/$sample.variants.chr${chr}.raw.all.vcf $output/$sample.variants.chr${chr}.raw.all.pbt.vcf $run_info
				fi

                ### filter this file to keep only the variants calls
                cat $output/$sample.variants.chr${chr}.raw.all.vcf | awk '$5 != "." || $0 ~ /^#/' | grep -v "\./\."  | grep -v "0\/0" > $output/$sample.variants.chr${chr}.raw.vcf
                sed '/^$/d' $output/$sample.variants.chr${chr}.raw.vcf > $output/$sample.variants.chr${chr}.raw.vcf.temp
                mv $output/$sample.variants.chr${chr}.raw.vcf.temp $output/$sample.variants.chr${chr}.raw.vcf

                ### prepare the file for backfilling
                cat $output/$sample.variants.chr${chr}.raw.all.vcf | grep -v "\./\." > $output/$sample.variants.chr${chr}.raw.all.vcf.temp
                rm $output/$sample.variants.chr${chr}.raw.all.vcf.idx
                mv $output/$sample.variants.chr${chr}.raw.all.vcf.temp $output/$sample.variants.chr${chr}.raw.all.vcf
                sed '/^$/d' $output/$sample.variants.chr${chr}.raw.all.vcf > $output/$sample.variants.chr${chr}.raw.all.vcf.temp
                mv $output/$sample.variants.chr${chr}.raw.all.vcf.temp $output/$sample.variants.chr${chr}.raw.all.vcf
                $tabix/bgzip $output/$sample.variants.chr${chr}.raw.all.vcf
            else
                param="-L chr${chr}"
                bam="-I $output/$sample.chr${chr}-sorted.bam"
                $script_path/unifiedgenotyper.sh "$bam" $output/$sample.variants.chr${chr}.raw.vcf BOTH "$param" EMIT_VARIANTS_ONLY $run_info
                rm $output/$sample.variants.chr${chr}.raw.vcf.idx

				## add phase by transmission if pedigree information provided.
				if [ $ped != "NA" ]
				then
					$script_path/phaseByTransmission.sh $output/$sample.variants.chr${chr}.raw.all.vcf $output/$sample.variants.chr${chr}.raw.all.pbt.vcf $run_info
				fi
            fi
        elif [ $SNV_caller == "SNVMIX" ]
		then
            ### call indels using GATK
            if [[ $all_sites == "YES"  && $tool == "exome" ]]
            then
                cat $TargetKit | grep -w chr$chr > $output/$sample.$chr.target.bed
                len=`cat $output/$sample.$chr.target.bed |wc -l`
                if [[ $only_ontarget == "YES" && $len -gt 0 ]]
                then
                    param="-L $output/$sample.$chr.target.bed"
                else
                    param="-L chr${chr}"
                fi
                ## call indels
                bam="-I $output/$sample.chr${chr}-sorted.bam"
                $script_path/unifiedgenotyper.sh "$bam" $output/$sample.variants.chr${chr}.raw.indel.all.vcf INDEL "$param" EMIT_ALL_SITES $run_info
                ### call snvs using SNVmix
                $script_path/snvmix2.sh $sample $input/chr${chr}.pileup $output/$sample.variants.chr${chr}.raw.snv.all.vcf all "$param" $run_info
                ### annoatte vcf
                $script_path/annotate_vcf.sh $output/$sample.variants.chr${chr}.raw.snv.all.vcf $chr $run_info "$bam"
                ### merge snvs and indels to give on vcf
                in="-V $output/$sample.variants.chr${chr}.raw.snv.all.vcf -V $output/$sample.variants.chr${chr}.raw.indel.all.vcf "
                $script_path/combinevcf.sh "$in" $output/$sample.variants.chr${chr}.raw.all.vcf $run_info yes
				in="-V $output/$sample.variants.chr${chr}.raw.snv.vcf.multi.vcf -V $output/$sample.variants.chr${chr}.raw.indel.vcf.multi.vcf"
                $script_path/combinevcf.sh "$in" $output/$sample.variants.chr${chr}.raw.multi.vcf $run_info yes
                cat $output/$sample.variants.chr${chr}.raw.all.vcf | awk '$5 != "." || $0 ~ /^#/' | grep -v "\./\."  | grep -v "0\/0" > $output/$sample.variants.chr${chr}.raw.vcf
                rm $output/$sample.variants.chr${chr}.raw.all.vcf.idx
                $tabix/bgzip $output/$sample.variants.chr${chr}.raw.all.vcf
				rm $output/$sample.$chr.target.bed
			else
                ## call indeles using GATK
                param="-L chr${chr}"
                bam="-I $output/$sample.chr${chr}-sorted.bam"
                $script_path/unifiedgenotyper.sh "$bam" $output/$sample.variants.chr${chr}.raw.indel.vcf INDEL "$param" EMIT_VARIANTS_ONLY $run_info
                ### call snvs using snvmix
                $script_path/snvmix2.sh $sample $input/chr${chr}.pileup $output/$sample.variants.chr${chr}.raw.snv.vcf target "$param" $run_info
                ## annotaet vcf
                $script_path/annotate_vcf.sh $output/$sample.variants.chr${chr}.raw.snv.vcf $chr $run_info "$bam"
                ### merge snvs and indels to give on vcf
                in="-V $output/$sample.variants.chr${chr}.raw.snv.vcf -V $output/$sample.variants.chr${chr}.raw.indel.vcf"
                $script_path/combinevcf.sh "$in" $output/$sample.variants.chr${chr}.raw.vcf $run_info yes
				in="-V $output/$sample.variants.chr${chr}.raw.snv.vcf.multi.vcf -V $output/$sample.variants.chr${chr}.raw.indel.vcf.multi.vcf"
                $script_path/combinevcf.sh "$in" $output/$sample.variants.chr${chr}.raw.multi.vcf $run_info yes
			fi
		fi
	else
        ## assuming that normal is the first column/sample
        normal=${sampleArray[1]}.chr$chr-sorted.bam
        inputfiles="-I $output/$normal ";

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
            fi
	    ### annotate vcfs
            in="-I $input/$bam"
            $script_path/annotate_vcf.sh $output/$sample.chr$chr.snv.vcf $chr $run_info "$in"
	    ### somatic indel calling
            $script_path/somaticindel.sh $output/$tumor $output/$normal $chr $sample $output $sample.chr$chr.indel.vcf $run_info
            ### annoatte vcfs
            $script_path/annotate_vcf.sh $output/$sample.chr$chr.indel.vcf $chr $run_info "$in"
        done
        param="-L chr${chr}"
		$script_path/unifiedgenotyper.sh "$inputfiles" $output/variants.chr${chr}.raw.vcf BOTH "$param" EMIT_VARIANTS_ONLY $run_info
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
        in="-V $output/MergeAllSamples.chr$chr.snvs.raw.vcf -V $output/MergeAllSamples.chr$chr.Indels.raw.vcf"
        $script_path/combinevcf.sh "$in" $output/MergeAllSamples.chr$chr.raw.vcf $run_info yes
		in="-V $output/MergeAllSamples.chr$chr.snvs.raw.multi.vcf -V $output/MergeAllSamples.chr$chr.Indels.raw.multi.vcf"
		$script_path/combinevcf.sh "$in" $output/MergeAllSamples.chr$chr.raw.multi.vcf $run_info yes
    fi

    ## remove files
    if [ ${#sampleArray[@]} == 1 ]
    then
        bam=chr${chr}.cleaned.bam
        rm $output/$bam.$chr.bam
        rm $output/$bam.$chr.bam.bai
        rm $output/${sampleArray[1]}.chr$chr.bam
        rm $output/${sampleArray[1]}.chr$chr.bam.bai
        rm $output/${sampleArray[1]}.chr$chr-sorted.bam
        rm $output/${sampleArray[1]}.chr$chr-sorted.bam.bai
    fi
    if [ ${#sampleArray[@]} -gt 1 ]
    then
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
    fi

    ## update dash board
    if [ $SGE_TASK_ID == 1 ]
    then
        $script_path/dashboard.sh $samples $run_info VariantCalling complete
    fi
    echo `date`
fi
