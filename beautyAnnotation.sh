#!/bin/sh
USR='TU03318'
WRDPS='bsu96681'

if [ $# != 4 ];
then
    echo "Usage <input .vcf file> <gene reference .txt file> <output file name>";
else
	set -x
	echo `date`
	misc=$1
    input=$2
    sample=$3
    run_info=$4
	
	tool_info=$( cat $run_info | grep -w '^TOOL_INFO' | cut -d '=' -f2)
	java=$( cat $tool_info | grep -w '^JAVA' | cut -d '=' -f2)
	REF=$( cat $tool_info | grep -w '^UCSC_REF_FLAT' | cut -d '=' -f2)
	chr=$(cat $run_info | grep -w '^CHRINDEX' | cut -d '=' -f2 | tr ":" "\n" | head -n $SGE_TASK_ID | tail -n 1)
	
	### SNVs
	
	INPUT=$input/$sample.variants.chr$chr.SNV.filter.i.c.vcf
	OUTPUT=$misc/$sample.variants.chr$chr.SNV.filter.i.c.misc.vcf
	$java/java -Xmx2g -Xms512m -jar $script_path/gene_annot.jar $REF $INPUT;
	
	NEWIN=${INPUT%.*}
	NEWIN=$NEWIN"_gene.annot";
	## additional parameters to prevent header (--header 0)
	if [ ! -f $NEWIN ]
	then
		echo ""
		echo "ERROR: Unable to map to genes, please contact Raymond Moore <moore.raymond@mayo.edu>"
		echo ""
		exit 1;
	fi
	$java/java -Xmx2g -Xms512m -jar $script_path/DrugTragetability.jar -u $USR -p $WRDPS -i $NEWIN -o $OUTPUT;
	## rm $NEWIN
	### INDELS
	INPUT=$input/$sample.variants.chr$chr.INDEL.filter.i.c.vcf
	OUTPUT=$input/$sample.variants.chr$chr.INDEL.filter.i.c.misc.vcf
	$java/java -Xmx2g -Xms512m -jar $script_path/gene_annot.jar $REF $INPUT;
	
	NEWIN=${INPUT%.*}
	NEWIN=$NEWIN"_gene.annot";
	## additional parameters to prevent header (--header 0)
	if [ ! -f $NEWIN ]
	then
		echo ""
		echo "ERROR: Unable to map to genes, please contact Raymond Moore <moore.raymond@mayo.edu>"
		echo ""
		exit 1;
	fi
	$java/java -Xmx2g -Xms512m -jar $script_path/DrugTragetability.jar -u $USR -p $WRDPS -i $NEWIN -o $OUTPUT;
	echo `date`
fi