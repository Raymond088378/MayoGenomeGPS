#!/bin/bash

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

if [ $# -le 4 ];
then
    echo -e "script to run SIFT annotation tool on a vcf file\nUsage: ./sift.sh <sift dir> <input dir><sample><run info><somatic/germline><SGE_TASK_ID(optional)>";
else
    set -x
    echo `date` 
    sift=$1
    input=$2
    sample=$3
    run_info=$4
    type=`echo $5 | tr "[A-Z]" "[a-z]"`	
	if [ $type == "somatic" ]
	then
		prefix="TUMOR"
		sam=$prefix.$sample
	else
		sam=$sample
	fi	
	if [ $6 ]
    then
        SGE_TASK_ID=$6
    fi	
	tool_info=$( cat $run_info | grep -w '^TOOL_INFO' | cut -d '=' -f2)
    sift_ref=$( cat $tool_info | grep -w '^SIFT_REF' | cut -d '=' -f2) 
    sift_path=$( cat $tool_info | grep -w '^SIFT' | cut -d '=' -f2) 
    script_path=$( cat $tool_info | grep -w '^WORKFLOW_PATH' | cut -d '=' -f2 )
    chr=$(cat $run_info | grep -w '^CHRINDEX' | cut -d '=' -f2 | tr ":" "\n" | head -n $SGE_TASK_ID | tail -n 1)
    multi=$( cat $run_info | grep -w '^MULTISAMPLE' | cut -d '=' -f2| tr "[a-z]" "[A-Z]")
	analysis=$( cat $run_info | grep -w '^ANALYSIS' | cut -d '=' -f2 | tr "[A-Z]" "[a-z]" )
	if [ $analysis == "annotation" ]
	then
		previous="reformat_VARIANTs.sh"
	else
		previous="OnTarget_variant.sh"
	fi	
	### update dashboard
	if [ $multi == "YES" ]
	then
		sample_info=$( cat $run_info | grep -w '^SAMPLE_INFO' | cut -d '=' -f2)
		ss=$( cat $sample_info | grep -w '^$sample' | cut -d '=' -f2 | tr "\t" " ")
		for i in $ss
		do
			$script_path/dashboard.sh $i $run_info Annotation started
		done
	else
		$script_path/dashboard.sh $sample $run_info Annotation started
	fi	
    ### hard coded
    snv_file=$sam.variants.chr$chr.SNV.final.i.c.vcf
    if [ ! -s $input/$snv_file ]
    then
		touch $input/$snv_file.sift.fix.log
        $script_path/email.sh $input/$snv_file "not found" $previous $run_info
		$script_path/wait.sh $input/$snv_file.sift.fix.log
    fi
    
    num_snvs=`cat $input/$snv_file | awk '$0 !~ /^#/' | wc -l`
	$script_path/filesize.sh sift $sam $input $snv_file $run_info
    #sift acceptable format 
    
    if [[ $num_snvs == 0 || $chr == 'M' ]]
    then
        touch $sift/${sam}_chr${chr}_predictions.tsv
        echo -e "Coordinates\tCodons\tTranscript ID\tProtein ID\tSubstitution\tRegion\tdbSNP ID\tSNP Type\tPrediction\tScore\tMedian Info\t# Seqs at position\tGene ID\tGene Name\tOMIM Disease\tAverage Allele Freqs\tUser Comment" > $sift/${sam}_chr${chr}_predictions.tsv
    else
        cat $input/$snv_file | awk '$0 !~ /^#/' | awk '$5 !~ /,/' | sed -e '/chr/s///g' | awk '{print $1","$2",1,"$4"/"$5}' > $sift/$snv_file.sift
        a=`pwd`
        #running SIFT for each sample
        cd $sift_path
        perl $sift_path/SIFT_exome_nssnvs.pl -i $sift/$snv_file.sift -d $sift_ref -o $sift/ -A 1 -B 1 -J 1 -K 1 > $sift/$sam.chr${chr}.sift.run
        id=`perl -n -e ' /Your job id is (\d+)/ && print "$1\n" ' $sift/$sam.chr${chr}.sift.run`
        rm $sift/$sam.chr${chr}.sift.run
        # sift inconsistent results flips alt base by itself getting rid of wrong calls from sift output
        mv $sift/$id/${id}_predictions.tsv $sift/${sam}_chr${chr}_predictions.tsv
        if [ ${#id} -gt 1 ]
		then
			rm -R $sift/$id
        fi
		$script_path/sift.inconsistent.pl $sift/${sam}_chr${chr}_predictions.tsv $sift/$snv_file.sift
        mv $sift/${sam}_chr${chr}_predictions.tsv_mod $sift/${sam}_chr${chr}_predictions.tsv
        rm $sift/$snv_file.sift
        cd $a
    fi
    if [ ! -s $sift/${sam}_chr${chr}_predictions.tsv ]
    then
        $script_path/errorlog.sh $sift/${sam}_chr${chr}_predictions.tsv sift.sh ERROR "failed to create" 
		exit 1;
    else
		$script_path/filesize.sh sift.out $sam $sift ${sam}_chr${chr}_predictions.tsv $run_info
	fi
    echo `date`
fi	
