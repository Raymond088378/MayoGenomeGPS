

if [ $# != 5 ]
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
    #SGE_TASK_ID=2
    tool_info=$( cat $run_info | grep -w '^TOOL_INFO' | cut -d '=' -f2)
    tool=$( cat $run_info | grep -w '^TYPE' | cut -d '=' -f2|tr "[A-Z]" "[a-z]")
	tabix=$( cat $tool_info | grep -w '^TABIX' | cut -d '=' -f2)
	perllib=$( cat $tool_info | grep -w '^PERLLIB_VCF' | cut -d '=' -f2)
	samtools=$( cat $tool_info | grep -w '^SAMTOOLS' | cut -d '=' -f2)	
    ref=$( cat $tool_info | grep -w '^REF_GENOME' | cut -d '=' -f2)
    gatk=$( cat $tool_info | grep -w '^GATK' | cut -d '=' -f2)
    dbSNP=$( cat $tool_info | grep -w '^dbSNP_REF' | cut -d '=' -f2)
    Kgenome=$( cat $tool_info | grep -w '^KGENOME_REF' | cut -d '=' -f2)
    java=$( cat $tool_info | grep -w '^JAVA' | cut -d '=' -f2)
    picard=$( cat $tool_info | grep -w '^PICARD' | cut -d '=' -f2 ) 
    script_path=$( cat $tool_info | grep -w '^WHOLEGENOME_PATH' | cut -d '=' -f2 )
    somatic_sniper=$( cat $tool_info | grep -w '^SOMATIC_SNIPER' | cut -d '=' -f2 )
    snvmix=$( cat $tool_info | grep -w '^SNVmix' | cut -d '=' -f2)
	version=$( cat $run_info | grep -w '^VERSION' | cut -d '=' -f2)
	samtools=$( cat $tool_info | grep -w '^SAMTOOLS' | cut -d '=' -f2)
	analysis=$( cat $run_info | grep -w '^ANALYSIS' | cut -d '=' -f2 |tr "[A-Z]" "[a-z]")
	run_num=$( cat $run_info | grep -w '^OUTPUT_FOLDER' | cut -d '=' -f2)
    flowcell=`echo $run_num | awk -F'_' '{print $NF}' | sed 's/.\(.*\)/\1/'`
	dbSNP_rsIDs=$( cat $tool_info | grep -w '^dbSNP_SNV_rsIDs' | cut -d '=' -f2 )
    sampleNames=$( echo $samples | tr ":" "\n" )
	chr=$( cat $run_info | grep -w '^CHRINDEX' | cut -d '=' -f2 | tr ":" "\n" | head -n $SGE_TASK_ID | tail -n 1 )
    out=$( cat $run_info | grep -w '^BASE_OUTPUT_DIR' | cut -d '=' -f2)
    PI=$( cat $run_info | grep -w '^PI' | cut -d '=' -f2)
    tool=$( cat $run_info | grep -w '^TYPE' | cut -d '=' -f2 | tr "[A-Z]" "[a-z]" )    
	all_sites=$( cat $tool_info | grep -w '^EMIT_ALL_SITES' | cut -d '=' -f2 | tr "[a-z]" "[A-Z]" )
    depth_filter=$( cat $tool_info | grep -w '^DEPTH_FILTER' | cut -d '=' -f2)
    prob_filter=$( cat $tool_info | grep -w '^PROB_FILTER' | cut -d '=' -f2)
	only_ontarget=$( cat $tool_info | grep -w '^TARGETTED' | cut -d '=' -f2 | tr "[a-z]" "[A-Z]" )
    TargetKit=$( cat $tool_info | grep -w '^ONTARGET' | cut -d '=' -f2 )
	threads=$( cat $tool_info | grep -w '^THREADS' | cut -d '=' -f2)
	version=$( cat $run_info | grep -w '^VERSION' | cut -d '=' -f2)
	SNV_caller=$( cat $run_info | grep -w '^SNV_CALLER' | cut -d '=' -f2)
    bedtools=$( cat $tool_info | grep -w '^BEDTOOLS' | cut -d '=' -f2 )
    PATH=$bedtools/:$PATH
	export PERL5LIB=$perllib
	PATH=$tabix/:$PATH
	bam=chr${chr}.cleaned.bam
	if [ $SGE_TASK_ID == 1 ]
    then
        if [ $analysis == "mayo" -o $analysis == "realign-mayo" ]
        then
            s=`echo $samples | tr ":" " "`
            for sam in $s
            do
                pos=$( cat $run_info | grep -w '^SAMPLENAMES' | cut -d '=' -f2 | tr ":" "\n" | grep -n $sam | cut -d ":" -f1)
                lanes=$( cat $run_info | grep -w '^LANEINDEX' | cut -d '=' -f2 | tr ":" "\n" | head -n $pos | tail -n 1 | tr "," " ")
				i=1
				for lane in $lanes
				do
					index=$( cat $run_info | grep -w '^LABINDEXES' | cut -d '=' -f2 | tr ":" "\n" | head -n $pos | tail -n 1 | tr "," "\n" | head -n $i | tail -n 1)
					if [ $index == "-" ]
					then
						$java/java -jar $script_path/AddSecondaryAnalysis.jar -p $script_path/AddSecondaryAnalysis.properties -l $lane -f $flowcell -r $run_num -s VariantCalling -a WholeGenome -v $version
					else
						$java/java -jar $script_path/AddSecondaryAnalysis.jar -p $script_path/AddSecondaryAnalysis.properties -l $lane -f $flowcell -i $index -r $run_num -s VariantCalling -a WholeGenome -v $version
					fi
					let i=i+1
				done		
            done
        fi
    fi    

	i=1
   for sample in $sampleNames
   do
       sampleArray[$i]=$sample
       let i=i+1
   done
    
    num_samples=`echo $samples | tr ":" " "`
#   gr=""
#   for i in $num_sample
#   do
#       a="grep -v \"^@RG.ID:$i\"";
#       gr="$gr| $a"
#    done    
    ## if multiple samples then split using read group and Validate BAM other wise just check and validate BAM

    if [ ! -s $input/$bam ]
    then
        echo "ERROR : variants.sh File $input/$bam does not exist"
        exit 1
    fi
    
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
            $samtools/samtools view -H $output/$sample.chr$chr.rg.bam | grep -E -v '$gr' | $samtools/samtools reheader - $output/$sample.chr$chr.rg.bam > $output/$sample.chr$chr.rg.re.bam
            mv $output/$sample.chr$chr.rg.re.bam $output/$sample.chr$chr.rg.bam

            if [ ! -s $output/$sample.chr$chr.rg.bam ]
            then
                echo "ERROR : variant_calls_per_chr File $output/$sample.chr$chr.rg.bam was not generated"
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
        
				$java/java -Xmx3g -Xms512m -jar $gatk/GenomeAnalysisTK.jar \
				-R $ref \
				-et NO_ET \
				-T UnifiedGenotyper \
				-glm BOTH \
				-nt $threads \
				$param \
				--output_mode EMIT_ALL_SITES \
				-I $output/$sample.chr${chr}-sorted.bam \
				--dbsnp $dbSNP \
				--out $output/$sample.variants.chr${chr}.raw.all.vcf
                rm $output/$sample.$chr.target.bed

				### filter this file to keep only the variants calls
                ### generate the raw calls
                cat $output/$sample.variants.chr${chr}.raw.all.vcf | awk '$5 != "." || $0 ~ /^#/' | grep -v "\./\." > $output/$sample.variants.chr${chr}.raw.vcf 
                sed '/^$/d' $output/$sample.variants.chr${chr}.raw.vcf > $output/$sample.variants.chr${chr}.raw.vcf.temp
                mv $output/$sample.variants.chr${chr}.raw.vcf.temp $output/$sample.variants.chr${chr}.raw.vcf
                
                ### prepare the file for backfilling
                cat $output/$sample.variants.chr${chr}.raw.all.vcf | grep -v "\./\." > $output/$sample.variants.chr${chr}.raw.all.vcf.temp
                rm $output/$sample.variants.chr${chr}.raw.all.vcf.idx
                mv $output/$sample.variants.chr${chr}.raw.all.vcf.temp $output/$sample.variants.chr${chr}.raw.all.vcf
                sed '/^$/d' $output/$sample.variants.chr${chr}.raw.all.vcf > $output/$sample.variants.chr${chr}.raw.all.vcf.temp
                mv $output/$sample.variants.chr${chr}.raw.all.vcf.temp $output/$sample.variants.chr${chr}.raw.all.vcf
				if [ $only_ontarget == "YES" ]
                then
                    perl $script_path/filter.dp.vcf.pl -i $output/$sample.variants.chr${chr}.raw.all.vcf -s $sample -d $depth_filter -o $output/$sample.variants.chr${chr}.raw.all.vcf.filter
                    mv $output/$sample.variants.chr${chr}.raw.all.vcf.filter $output/$sample.variants.chr${chr}.raw.all.vcf 
                fi
				sed '/^$/d' $output/$sample.variants.chr${chr}.raw.all.vcf > $output/$sample.variants.chr${chr}.raw.all.vcf.temp
				mv $output/$sample.variants.chr${chr}.raw.all.vcf.temp $output/$sample.variants.chr${chr}.raw.all.vcf
                $tabix/bgzip $output/$sample.variants.chr${chr}.raw.all.vcf   
			else
				$java/java -Xmx3g -Xms512m -jar $gatk/GenomeAnalysisTK.jar \
				-R $ref \
				-et NO_ET \
				-T UnifiedGenotyper \
				-glm BOTH \
				-nt $threads \
				-L chr${chr} \
				-I $output/$sample.chr${chr}-sorted.bam \
				--dbsnp $dbSNP \
				--out $output/$sample.variants.chr${chr}.raw.vcf 
                rm $output/$sample.variants.chr${chr}.raw.vcf.idx 
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

                $java/java -Xmx3g -Xms512m -jar $gatk/GenomeAnalysisTK.jar \
                -R $ref \
                -et NO_ET \
                -T UnifiedGenotyper \
                --output_mode EMIT_ALL_SITES \
                -glm INDEL \
                -nt $threads \
                $param \
                -I $output/$sample.chr${chr}-sorted.bam \
                --dbsnp $dbSNP \
                --out $output/$sample.variants.chr${chr}.raw.indel.all.vcf 

				### call snvs using SNVmix
				$snvmix/SNVMix2 -i $input/chr${chr}.pileup -f -m $snvmix/Mu_pi.txt -o $output/$sample.variants.chr${chr}.raw.snv.all
				if [ $only_ontarget == "YES" ]
                then
                    perl $script_path/snvmix_to_vcf.pl -i $output/$sample.variants.chr${chr}.raw.snv.all -o $output/$sample.variants.chr${chr}.raw.snv.all.vcf -s $sample -p $prob_filter -d $depth_filter
                    len=`cat $output/$sample.$chr.target.bed |wc -l`
					if [ $len -gt 0 ]
					then
						$bedtools/intersectBed -a $output/$sample.variants.chr${chr}.raw.snv.all.vcf -b $output/$sample.$chr.target.bed -wa -header > $output/$sample.variants.chr${chr}.raw.snv.all.vcf.i
                    else
						cp $output/$sample.variants.chr${chr}.raw.snv.all.vcf $output/$sample.variants.chr${chr}.raw.snv.all.vcf.i	
					fi	
					mv $output/$sample.variants.chr${chr}.raw.snv.all.vcf.i $output/$sample.variants.chr${chr}.raw.snv.all.vcf
                    rm $output/$sample.$chr.target.bed
                else
                    perl $script_path/snvmix_to_vcf.pl -i $output/$sample.variants.chr${chr}.raw.snv.all -o $output/$sample.variants.chr${chr}.raw.snv.all.vcf -s $sample
                    rm $output/$sample.$chr.target.bed
                fi
				### annotate vcf from snvmix
				$java/java -Xmx3g -Xms512m -jar $gatk/GenomeAnalysisTK.jar \
				-R $ref \
				-et NO_ET \
				-T VariantAnnotator \
				-I $output/$sample.chr${chr}-sorted.bam \
				-V $output/$sample.variants.chr${chr}.raw.snv.all.vcf \
				--dbsnp $dbSNP \
				-L chr$chr \
				-A QualByDepth -A MappingQualityRankSumTest -A HaplotypeScore -A DepthOfCoverage -A MappingQualityZero -A RMSMappingQuality -A FisherStrand \
				--out $output/$sample.variants.chr${chr}.raw.snv.all.vcf.temp
				mv $output/$sample.variants.chr${chr}.raw.snv.all.vcf.temp $output/$sample.variants.chr${chr}.raw.snv.all.vcf
				mv $output/$sample.variants.chr${chr}.raw.snv.all.vcf.temp.idx $output/$sample.variants.chr${chr}.raw.snv.all.vcf.idx
				
				rm $output/$sample.variants.chr${chr}.raw.snv.all
				### merge snvs and indels to give on vcf
				$java/java -Xmx2g -Xms512m -jar $gatk/GenomeAnalysisTK.jar \
				-R $ref \
				-et NO_ET \
				-T CombineVariants \
				-V $output/$sample.variants.chr${chr}.raw.snv.all.vcf -V $output/$sample.variants.chr${chr}.raw.indel.all.vcf  \
				-o $output/$sample.variants.chr${chr}.raw.all.vcf

				if [ -s $output/$sample.variants.chr${chr}.raw.all.vcf ]
				then
					rm $output/$sample.variants.chr${chr}.raw.snv.all.vcf $output/$sample.variants.chr${chr}.raw.indel.all.vcf
                    rm $output/$sample.variants.chr${chr}.raw.snv.all.vcf.idx $output/$sample.variants.chr${chr}.raw.indel.all.vcf.idx    
				else
					echo "ERROR: merging indels and snvs failed for $sample (variants.sh)"
					exit 1;
				fi		
				cat $output/$sample.variants.chr${chr}.raw.all.vcf | awk '$5 != "N" || $0 ~ /^#/' | grep -v "\./\." > $output/$sample.variants.chr${chr}.raw.vcf
                rm $output/$sample.variants.chr${chr}.raw.all.vcf.idx
				$tabix/bgzip $output/$sample.variants.chr${chr}.raw.all.vcf	
			else
				## call indeles using GATK
				$java/java -Xmx3g -Xms512m -jar $gatk/GenomeAnalysisTK.jar \
				-R $ref \
				-et NO_ET \
				-T UnifiedGenotyper \
				-glm INDEL \
				-nt $threads \
				-L chr${chr} \
				-I $output/$sample.chr${chr}-sorted.bam \
				--dbsnp $dbSNP \
				--out $output/$sample.variants.chr${chr}.raw.indel.vcf
				
				### call snvs using snvmix
				$snvmix/SNVMix2 -i $input/chr${chr}.pileup -m $snvmix/Mu_pi.txt -o $output/$sample.variants.chr${chr}.raw.snv
				perl $script_path/snvmix_to_vcf.pl -i $output/$sample.variants.chr${chr}.raw.snv -o $output/$sample.variants.chr${chr}.raw.snv.vcf -s $sample
				### annotate vcf from snvmix
				$java/java -Xmx3g -Xms512m -jar $gatk/GenomeAnalysisTK.jar \
				-R $ref \
				-et NO_ET \
				-T VariantAnnotator \
				-I $output/$sample.chr${chr}-sorted.bam \
				-V $output/$sample.variants.chr${chr}.raw.snv.vcf \
				--dbsnp $dbSNP \
				-L chr$chr \
				-A QualByDepth -A MappingQualityRankSumTest -A HaplotypeScore -A DepthOfCoverage -A MappingQualityZero -A RMSMappingQuality -A FisherStrand \
				--out $output/$sample.variants.chr${chr}.raw.snv.vcf.temp
				mv $output/$sample.variants.chr${chr}.raw.snv.vcf.temp $output/$sample.variants.chr${chr}.raw.snv.vcf
				mv $output/$sample.variants.chr${chr}.raw.snv.vcf.temp.idx $output/$sample.variants.chr${chr}.raw.snv.vcf.idx	
				
				rm $output/$sample.variants.chr${chr}.raw.snv
				### merge snvs and indels to give on vcf
				$java/java -Xmx2g -Xms512m -jar $gatk/GenomeAnalysisTK.jar \
				-R $ref \
				-et NO_ET \
				-T CombineVariants \
				-V $output/$sample.variants.chr${chr}.raw.snv.vcf -V $output/$sample.variants.chr${chr}.raw.indel.vcf \
				-o $output/$sample.variants.chr${chr}.raw.vcf 
				
				if [ -s $output/$sample.variants.chr${chr}.raw.vcf ]
				then
					rm $output/$sample.variants.chr${chr}.raw.snv.vcf $output/$sample.variants.chr${chr}.raw.indel.vcf
                    rm $output/$sample.variants.chr${chr}.raw.snv.vcf.idx $output/$sample.variants.chr${chr}.raw.indel.vcf.idx
				else
					echo "ERROR: merging indels and snvs failed for $sample (variant_calls_per_chr.sh)"
				fi	
			fi		
		fi
		if [ ! -s $output/$sample.variants.chr${chr}.raw.vcf ]
        then
            echo "ERROR : variants.sh File $output/$sample.variants.chr${chr}.raw.vcf not generated "
            exit 1
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

            ##run somatic sniper 
            $somatic_sniper/bam-somaticsniper -q 20 -Q 30 -f $ref $output/$tumor $output/$normal $output/$snv     

            if [ ! -s $output/$snv ]
            then		
                echo "ERROR :variants.sh SomaticSnipper failed, file $output/$snv not generated "
                exit 1
            fi
            
            ## convert sniper output to VCF
            perl $script_path/ss2vcf.pl $output/$snv $output/$sample.chr$chr.snv.vcf $dbSNP_rsIDs ${sampleArray[1]} $sample $output/$sample.chr$chr.snv.triallele.out
            
            if [ -s $output/$sample.chr$chr.snv.vcf ]
            then
                rm $output/$snv
            else
                echo "ERROR: $output/$sample.chr$chr.snv.vcf not found"
            fi
            
            # ## annotate SNVs
            $java/java -Xmx3g -Xms512m -jar $gatk/GenomeAnalysisTK.jar \
            -R $ref \
            -et NO_ET \
            -T VariantAnnotator \
            -I $output/$tumor \
            -I $output/$normal \
            -V $output/$sample.chr$chr.snv.vcf \
            --dbsnp $dbSNP \
            -L chr$chr \
            -A QualByDepth -A MappingQualityRankSumTest -A HaplotypeScore -A DepthOfCoverage -A MappingQualityZero -A RMSMappingQuality -A FisherStrand \
            --out $output/$sample.chr$chr.snv.vcf.temp	

            if [ -s $output/$sample.chr$chr.snv.vcf.temp ]
            then
                mv $output/$sample.chr$chr.snv.vcf.temp $output/$sample.chr$chr.snv.vcf
                rm $output/$sample.chr$chr.snv.triallele.out
            else		
                echo "ERROR : variants.sh SNV VariantAnnotator failed, file:$output/$sample.chr$chr.snv.vcf.temp"
                exit 1
            fi

            if [ -s $output/$sample.chr$chr.snv.vcf.temp.idx ]
            then
                mv $output/$sample.chr$chr.snv.vcf.temp.idx $output/$sample.chr$chr.snv.vcf.idx
            else	
                echo "ERROR: variants.sh SNV VariantAnnotator did not generated index, file: $output/$sample.chr$chr.snv.vcf.idx" 
                exit 1
            fi
            
            ## Somatic Indel detector
            indel=$sample.chr$chr.indel.vcf
            indel_v=$sample.chr$chr.indel.txt
            $java/java -Xmx3g -Xms512m -jar $gatk/GenomeAnalysisTK.jar \
            -R $ref \
            -et NO_ET \
            -T SomaticIndelDetector \
            --window_size 1000 \
            -o $output/$indel \
            -verbose $output/$indel_v \
            -I:normal $output/$normal \
            -I:tumor $output/$tumor

            if [ ! -s $output/$indel ]
            then
                echo "ERROR : variants.sh SomaticIndelDetector failed, file $output/$indel not generated "
                exit 1
            else
                rm $output/$indel_v
            fi
            
            # ## annotate INDELs
            $java/java -Xmx3g -Xms512m -jar $gatk/GenomeAnalysisTK.jar \
            -R $ref \
            -et NO_ET \
            -T VariantAnnotator \
            -I $output/$tumor \
            -I $output/$normal \
            -V $output/$indel \
            --dbsnp $dbSNP \
            -L chr$chr \
            -A QualByDepth -A MappingQualityRankSumTest -A HaplotypeScore -A DepthOfCoverage -A MappingQualityZero -A RMSMappingQuality -A FisherStrand \
            --out $output/$indel.temp

            if [ -s $output/$indel.temp ]
            then
                mv $output/$indel.temp $output/$indel
            else
                echo "ERROR : variants.sh Indel VariantAnnotator failed, file:$output/$indel.temp"
                exit 1;
            fi

            if [ -s $output/$indel.temp.idx ]
            then
                mv $output/$indel.temp.idx $output/$indel.idx	
            else	
                echo "ERROR: variants.sh SNV VariantAnnotator did not generated index, file: $output/$indel.temp.idx"
                exit 1;
            fi
        done	

    #    ##run GATK unified genotyper on all samples of a group
        $java/java -Xmx3g -Xms512m -jar $gatk/GenomeAnalysisTK.jar \
        -R $ref \
        -et NO_ET \
        -T UnifiedGenotyper \
		-nt $threads \
        -glm BOTH \
        -L chr${chr} \
        $inputfiles \
        --out $output/variants.chr${chr}.raw.vcf 
    
        if [ ! -s $output/variants.chr${chr}.raw.vcf ]
        then		
            echo "ERROR : variants.sh Unified Genotyper failed, file $output/variants.chr${chr}.raw.vcf not generated "
        fi
    fi
    ### after this we get multiple indels and snp files which need to be merged for Multi samples but just filter for 	
    if [ ${#sampleArray[@]} -gt 1 ]
    then
        ### INDEL
        input_var=""
        input_files=""
        for i in $(seq 2 ${#sampleArray[@]})
        do
            sample=${sampleArray[$i]}
            indel=$sample.chr$chr.indel.vcf
            input_var="${input_var} -V $output/$indel"
            input_files="${input_files} $output/$indel"
            input_indexes="${input_indexes} $output/$indel.idx"
        done

        $java/java -Xmx3g -Xms512m -jar $gatk/GenomeAnalysisTK.jar \
        -R $ref \
        -et NO_ET \
        -T CombineVariants \
        $input_var \
        -o  $output/MergeAllSamples.chr$chr.Indels.raw.vcf	

        if [ ! -s $output/MergeAllSamples.chr$chr.Indels.raw.vcf ]
        then		
            echo "ERROR : variants.sh CombineVariants indels failed, file $output/MergeAllSamples.chr$chr.Indels.raw.vcf not generated "
            exit 1;
        else
            rm $input_files
            rm $input_indexes
        fi

        ##Merge SNVs 
        input_var=""
        input_files=""
        for i in $(seq 2 ${#sampleArray[@]})
        do
            sample=${sampleArray[$i]}
            snv=$sample.chr$chr.snv.vcf
            input_var="${input_var} -V $output/$snv"
            input_files="${input_files} $output/$snv"
            input_indexes="${input_indexes} $output/$snv.idx"
        done

        $java/java -Xmx3g -Xms512m -jar $gatk/GenomeAnalysisTK.jar \
        -R $ref \
        -et NO_ET \
        -T CombineVariants \
        $input_var \
        -o  $output/MergeAllSamples.chr$chr.snvs.raw.vcf	

        if [ ! -s $output/MergeAllSamples.chr$chr.snvs.raw.vcf ]
        then	
            echo "ERROR : variants.sh CombineVariants snv failed, file $output/MergeAllSamples.chr$chr.snvs.raw.vcf not generated "
            exit 1;
        else
            rm $input_files
            rm $input_indexes
        fi

        $java/java -Xmx3g -Xms512m -jar $gatk/GenomeAnalysisTK.jar \
        -R $ref \
        -et NO_ET \
        -T CombineVariants \
        -V $output/MergeAllSamples.chr$chr.snvs.raw.vcf -V $output/MergeAllSamples.chr$chr.Indels.raw.vcf \
        -o $output/MergeAllSamples.chr$chr.raw.vcf	

        if [ ! -s $output/MergeAllSamples.chr$chr.raw.vcf ]
        then		
            echo "ERROR : variants.sh CombineVariants Indel, Snv failed, file $output/MergeAllSamples.chr$chr.raw.vcf not generated "
            exit 1;
        else
            rm $output/MergeAllSamples.chr$chr.snvs.raw.vcf $output/MergeAllSamples.chr$chr.Indels.raw.vcf
            rm $output/MergeAllSamples.chr$chr.snvs.raw.vcf.idx $output/MergeAllSamples.chr$chr.Indels.raw.vcf.idx
        fi    
    fi
            
    
    ## remove files
    if [ ${#sampleArray[@]} == 1 ]
    then
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
            rm  $output/${sampleArray[$i]}.chr$chr-sorted.bam
            rm  $output/${sampleArray[$i]}.chr$chr-sorted.bam.bai
        done	
    fi
	if [ $SGE_TASK_ID == 1 ]
    then
        if [ $analysis == "mayo" -o $analysis == "realign-mayo" ]
        then
            s=`echo $samples | tr ":" " "`
            for sam in $s
            do
                pos=$( cat $run_info | grep -w '^SAMPLENAMES' | cut -d '=' -f2 | tr ":" "\n" | grep -n $sample | cut -d ":" -f1)
                lanes=$( cat $run_info | grep -w '^LANEINDEX' | cut -d '=' -f2 | tr ":" "\n" | head -n $pos | tail -n 1 | tr "," " ")
				i=1
				for lane in $lanes
				do
					index=$( cat $run_info | grep -w '^LABINDEXES' | cut -d '=' -f2 | tr ":" "\n" | head -n $pos | tail -n 1 | tr "," "\n" | head -n $i | tail -n 1)
					if [ $index == "-" ]
					then
						$java/java -jar $script_path/AddSecondaryAnalysis.jar -p $script_path/AddSecondaryAnalysis.properties -l $lane -c -f $flowcell -r $run_num -s VariantCalling -a WholeGenome -v $version
					else
						$java/java -jar $script_path/AddSecondaryAnalysis.jar -p $script_path/AddSecondaryAnalysis.properties -l $lane -c -f $flowcell -i $index -r $run_num -s VariantCalling -a WholeGenome -v $version
					fi
					let i=i+1
				done			
            done
        fi
    fi    
	echo `date`
fi		
