#!/bin/sh

##	INFO
##	script is used to submit batch scripts to locally installed sift database

###############################
#	$1		=		sift output directory	
#	$2		=		sample name
#	$3		=		SNV input file
#	$4		=		directory for input file
#	$5		=		chromosome
#	$6		=		run_info file
#################################

if [ $# != 4 ];
then
    echo "Usage:<sift dir> <input dir><sample><run info>";
else
    set -x
    echo `date` 
    sift=$1
    input=$2
    sample=$3
    run_info=$4
    #SGE_TASK_ID=1
    tool_info=$( cat $run_info | grep -w '^TOOL_INFO' | cut -d '=' -f2)
    version=$( cat $run_info | grep -w '^VERSION' | cut -d '=' -f2)
    sift_ref=$( cat $tool_info | grep -w '^SIFT_REF' | cut -d '=' -f2) 
    sift_path=$( cat $tool_info | grep -w '^SIFT' | cut -d '=' -f2) 
    script_path=$( cat $tool_info | grep -w '^WHOLEGENOME_PATH' | cut -d '=' -f2 )
    analysis=$( cat $run_info | grep -w '^ANALYSIS' | cut -d '=' -f2| tr "[A-Z]" "[a-z]")
    run_num=$( cat $run_info | grep -w '^OUTPUT_FOLDER' | cut -d '=' -f2)
    flowcell=`echo $run_num | awk -F'_' '{print $NF}' | sed 's/.\(.*\)/\1/'`
    chr=$(cat $run_info | grep -w '^CHRINDEX' | cut -d '=' -f2 | tr ":" "\n" | head -n $SGE_TASK_ID | tail -n 1)
    java=$( cat $tool_info | grep -w '^JAVA' | cut -d '=' -f2 )
	
    ### update dashboard
    $script_path/dashboard.sh $sample $run_info Annotation started
    
    ### hard coded
    snv_file=$sample.variants.chr$chr.SNV.filter.i.c.vcf
	
    num_snvs=`cat $input/$snv_file | awk '$0 !~ /^#/' | wc -l`
    #sift acceptable format 
    
    if [ $num_snvs == 0 ]
    then
        touch $sift/${sample}_chr${chr}_predictions.tsv
        echo -e "Coordinates\tCodons\tTranscript ID\tProtein ID\tSubstitution\tRegion\tdbSNP ID\tSNP Type\tPrediction\tScore\tMedian Info\t# Seqs at position\tGene ID\tGene Name\tOMIM Disease\tAverage Allele Freqs\tUser Comment" > $sift/${sample}_chr${chr}_predictions.tsv
    else
        cat $input/$snv_file | awk '$0 !~ /^#/' | sed -e '/chr/s///g' | awk '{print $1","$2",1,"$4"/"$5}' > $sift/$snv_file.sift
        a=`pwd`
        #running SIFT for each sample
        cd $sift_path
        perl $sift_path/SIFT_exome_nssnvs.pl -i $sift/$snv_file.sift -d $sift_ref -o $sift/ -A 1 -B 1 -J 1 -K 1 > $sift/$sample.chr${chr}.sift.run
        id=`perl -n -e ' /Your job id is (\d+)/ && print "$1\n" ' $sift/$sample.chr${chr}.sift.run`
        rm $sift/$sample.chr${chr}.sift.run
        # sift inconsistent results flips alt base by itself getting rid of wrong calls from sift output
        mv $sift/$id/${id}_predictions.tsv $sift/${sample}_chr${chr}_predictions.tsv
        if [ ${#id} -gt 1 ]
		then
			rm -R $sift/$id
        fi
		perl $script_path/sift.inconsistent.pl $sift/${sample}_chr${chr}_predictions.tsv $sift/$snv_file.sift
        mv $sift/${sample}_chr${chr}_predictions.tsv_mod $sift/${sample}_chr${chr}_predictions.tsv
        rm $sift/$snv_file.sift
        cd $a
    fi
    if [ ! -s $sift/${sample}_chr${chr}_predictions.tsv ]
    then
        echo " ERROR : sift failed for $sample $chr "
    fi
    echo `date`
fi	
