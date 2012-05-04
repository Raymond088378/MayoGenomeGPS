#!/bin/sh

if [ $# != 4 ]
then
	echo "Usage: <raw vcf complete path><outputfile><type={BOTH,SNP,INDEL}><run info>"
else
    set -x
    echo `date`
    inputvcf=$1
    outfile=$2
    variant=$3
    run_info=$4
    
    input_dir=`dirname $inputvcf`
    input_name=`basename $inputvcf`
    output=`dirname $outfile`
	
	raw_calls=`cat $inputvcf | awk '$0 !~ /#/' | wc -l`
	
    if [ -d $output/temp ]
    then
            echo "temp already there"
    else
            mkdir $output/temp
    fi
    if [ -d $output/plot ]
    then
            echo "plot already there"
    else
            mkdir $output/plot
    fi

    tool_info=$( cat $run_info | grep -w '^TOOL_INFO' | cut -d '=' -f2)
    samtools=$( cat $tool_info | grep -w '^SAMTOOLS' | cut -d '=' -f2)	
    ref=$( cat $tool_info | grep -w '^REF_GENOME' | cut -d '=' -f2)
    gatk=$( cat $tool_info | grep -w '^GATK' | cut -d '=' -f2)
    dbSNP=$( cat $tool_info | grep -w '^dbSNP_REF' | cut -d '=' -f2)
    Kgenome=$( cat $tool_info | grep -w '^KGENOME_REF' | cut -d '=' -f2)
    hapmap=$( cat $tool_info | grep -w '^HAPMAP_VCF' | cut -d '=' -f2)
    omni=$( cat $tool_info | grep -w '^OMNI_VCF' | cut -d '=' -f2)
    java=$( cat $tool_info | grep -w '^JAVA' | cut -d '=' -f2)
    picard=$( cat $tool_info | grep -w '^PICARD' | cut -d '=' -f2 ) 
	script_path=$( cat $tool_info | grep -w '^WHOLEGENOME_PATH' | cut -d '=' -f2 )

    out=$( cat $run_info | grep -w '^BASE_OUTPUT_DIR' | cut -d '=' -f2)
    PI=$( cat $run_info | grep -w '^PI' | cut -d '=' -f2)
    tool=$( cat $run_info | grep -w '^TYPE' | cut -d '=' -f2 | tr "[A-Z]" "[a-z]" )
    run_num=$( cat $run_info | grep -w '^OUTPUT_FOLDER' | cut -d '=' -f2)


    if [ ! -s $hapmap ]
        then
        echo "ERROR: filter_variant_per_chr File $hapmap not available" 
        exit 1
    fi

    if [ ! -s $omni ]
        then
       
        echo "ERROR: filter_variant_vqsr File $omni not available" 
        exit 1
    fi

    if [ ! -s $dbSNP ]
        then
        echo "ERROR: filter_variant_vqsr File $dbsnp not available" 
        exit 1
    fi

    if [ ! -s $inputvcf ]
        then
        echo "ERROR: filter_variant_vqsr File $inputvcf not available"
        exit 1
    fi
    
    ## Variant Recalibrator SNP
    $java/java -Xmx3g -Xms512m -jar $gatk/GenomeAnalysisTK.jar \
    -R $ref \
    -et NO_ET \
    -T VariantRecalibrator \
    -mode $variant \
    -input $inputvcf \
    -resource:hapmap,known=false,training=true,truth=true,prior=15.0 $hapmap \
    -resource:omni,known=false,training=true,truth=false,prior=12.0 $omni \
    -resource:dbsnp,known=true,training=false,truth=false,prior=8.0 $dbSNP \
    -an QD -an HaplotypeScore -an MQRankSum -an ReadPosRankSum -an FS -an MQ \
    -recalFile $output/temp/$input_name.recal \
    -tranchesFile $output/temp/$input_name.tranches \
    --maxGaussians 4 \
    --percentBadVariants 0.05 \
    -rscriptFile $output/plot/$input_name.plots.R

    if [ ! -s $output/temp/$input_name.recal ]
        then
        echo "ERROR: filter_variant_per_chr, VariantRecalibrator File $output/temp/$input_name.recal not generated"
    fi

    if [ ! -s $output/temp/$input_name.tranches ]
        then
        echo "ERROR: filter_variant_per_chr, VariantRecalibrator File $output/temp/$input_name.tranches not generated"
    fi

    ## Apply Recalibrator
    $java/java -Xmx6g -XX:-UseGCOverheadLimit -Xms512m -jar $gatk/GenomeAnalysisTK.jar \
    -R $ref \
    -et NO_ET \
    -mode $variant \
    -T ApplyRecalibration \
    -input $inputvcf \
    -recalFile $output/temp/$input_name.recal \
    -tranchesFile $output/temp/$input_name.tranches \
    -o $outfile
	
	filter_calls=`cat $outfile | awk '$0 !~ /#/' | wc -l`
    if [[ ! -s $outfile || ! -s ${outfile}.idx || $raw_calls != $filter_calls ]]
    then
        echo "WARNING: filter_variant_vqsr, ApplyRecalibration File $outfile not generated"
        $java/java -Xmx1g -Xms512m -jar $gatk/GenomeAnalysisTK.jar \
        -R $ref \
        -et NO_ET \
        -l INFO \
        -T VariantFiltration \
        -V $inputvcf \
        -o $outfile \
        --clusterSize 10 \
        --filterExpression "MQ0 >= 4 && ((MQ0 / (1.0 * DP)) > 0.1)" \
        --filterName "HARD_TO_VALIDATE" \
        --filterExpression "QD < 1.5 " \
        --filterName "LowQD"
        if [ -s $outfile ]
        then
            echo "WARNING: VQSR failed so applied hard filters on the data"
        else
            echo "ERROR: hard filters failed to generate $outfile"
            exit 1;
        fi    
    fi
            
    ## remove the file
    if [ -s $output/temp/$input_name.recal ]
	then
		rm $output/temp/$input_name.recal
    fi
	
	if [ -s $output/temp/$input_name.tranches ]
	then
		rm $output/temp/$input_name.tranches	
    fi
	echo `date`
fi	