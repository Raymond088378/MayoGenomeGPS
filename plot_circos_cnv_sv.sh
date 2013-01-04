#!/bin/bash

if [ $# != 6 ]
then
    echo -e "script to make cirocs plot for structural variant and copy number variantion per sample\nUsage: <sv file break> <sv file crest><cnv file> <sample> <output dir> <run info>"
else
    set -x
    echo `date`
    sv_file_break=$1
    sv_file_crest=$2
    cnv_file=$3
    sample=$4
    out=$5
    run_info=$6
    
    tool_info=$( cat $run_info | grep -w '^TOOL_INFO' | cut -d '=' -f2)
    script_path=$( cat $tool_info | grep -w '^WORKFLOW_PATH' | cut -d '=' -f2)
    circos=$( cat $tool_info | grep -w '^CIRCOS' | cut -d '=' -f2)	
    perllib=$( cat $tool_info | grep -w '^PERLLIB_CIRCOS' | cut -d '=' -f2)	
	perl=$( cat $tool_info | grep -w '^PERL_CIRCOS' | cut -d '=' -f2)	
    multi=$( cat $run_info | grep -w '^MULTISAMPLE' | cut -d '=' -f2| tr "[a-z]" "[A-Z]")
	if [ $multi == "YES" ]
	then
		previous="summaryze_struct_group.sh"
	else	
		previous="summaryze_struct_single.sh"
	fi
    if [ ! -f $sv_file ]
    then
		$script_path/email.sh $sv_file "not exist" $previous $run_info
		touch $sv_file.fix.log
		$script_path/wait.sh $sv_file.fix.log
    fi    
    
    if [ ! -f $cnv_file ]
    then
		$script_path/email.sh $cnv_file "not exist" $previous $run_info
		touch $cnv_file.fix.log
		$script_path/wait.sh $cnv_file.fix.log
    fi
    
    cat $cnv_file | awk '$0 !~ /#/' | awk '{print $1"\t"$2"\t"$3"\t"$5}' | sort > $out/$sample.cnv.ForCircos.bed
    export PERL5LIB=$PERL5LIB:$perllib
    cat $sv_file_crest | awk '{print $1"\t"$2"\t"$3"\t"$5"\t"$6}' > $out/$sample.crest.bed
    cat $out/$sample.crest.bed $sv_file_break > $out/$sample_sv.bed
    rm $out/$sample.crest.bed
    $script_path/generate.sv_file.circos.pl $out/$sample_sv.bed > $out/$sample.sv.ForCircos.bed
    rm $out/$sample_sv.bed
    		
    ##script to create config files
    $script_path/circos_config.pl -p $script_path -v $out/$sample.sv.ForCircos.bed -c $out/$sample.cnv.ForCircos.bed -s $sample -o $out
    if [ ! -s $out/$sample.main.config ]
    then
        $script_path/errorlog.sh $out/$sample.main.config plot_circos_cnv_sv.sh ERROR "doesn't exist"
        exit 1;
    fi    
    ## plot circos
    $perl $circos/bin/circos -conf $out/$sample.main.config	
    ## deleting files
    rm $out/$sample.sv.ForCircos.bed
    rm  $out/$sample.cnv.ForCircos.bed
    rm $out/$sample.main.config
    echo `date`
fi
 
