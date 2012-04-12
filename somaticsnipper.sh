if [ $# != 8 ]
then
    echo "Usage: <normal bam> <tumor bam > <output dir> <chromosome> <tumor sample name> <normal sample name> <output vcf file name> <run info>"
else
    set -x
    echo `date`
    normal_bam=$1
    tumor_bam=$2
    output=$3
    chr=$4
    tumor_sample=$5
	normal_sample=$6
    output_file=$7
	run_info=$8
	
	tool_info=$( cat $run_info | grep -w '^TOOL_INFO' | cut -d '=' -f2)
	somatic_sniper=$( cat $tool_info | grep -w '^SOMATIC_SNIPER' | cut -d '=' -f2 )
	script_path=$( cat $tool_info | grep -w '^WHOLEGENOME_PATH' | cut -d '=' -f2 )
	ref=$( cat $tool_info | grep -w '^REF_GENOME' | cut -d '=' -f2)
	dbSNP_rsIDs=$( cat $tool_info | grep -w '^dbSNP_SNV_rsIDs' | cut -d '=' -f2 )
	snv=$tumor_sample.chr$chr.snv.output
	$somatic_sniper/bam-somaticsniper -q 20 -Q 30 -f $ref $tumor_bam $normal_bam $output/$snv   
	
	if [ ! -s $output/$snv ]
	then		
		echo "ERROR :variants.sh SomaticSnipper failed, file $output/$snv not generated "
		exit 1
	fi
	
	## convert sniper output to VCF
	perl $script_path/ss2vcf.pl $output/$snv $output/$output_file $dbSNP_rsIDs $normal_sample $tumor_sample $output/$output_file.triallele.out
	if [ -s $output/$output_file ]
	then
		rm $output/$snv
	else
		echo "ERROR: $output/$output_file not found"
		exit 1;
	fi
    echo `date`
fi	