#!/bin/sh

if [ $# != 2 ]
then
    echo "Usage : <output_dir> <run_info>";
else
    set -x
    echo `date`
    output_dir=$1
    run_info=$2
    tool_info=$( cat $run_info | grep -w '^TOOL_INFO' | cut -d '=' -f2)
    script_path=$( cat $tool_info | grep -w '^WHOLEGENOME_PATH' | cut -d '=' -f2 )
    java=$( cat $tool_info | grep -w '^JAVA' | cut -d '=' -f2 )
	genome_version=$(cat $run_info | grep -w '^GENOMEBUILD' | cut -d '=' -f2)
    ## jar script to add IGV, pathway, TIssue specificity and Gene Card link
    $java/java -Xmx2g -Xms512m -jar $script_path/exome_annot.jar annotate $tool_info $output_dir/Reports_per_Sample/ $genome_version
	
    cd $output_dir/Reports_per_Sample/
    ## sort the INDEL file
    for i in *.INDEL.*.xls
    do
        perl $script_path/sort.variantReport.pl -i $i -o $i.sort -f Start
        mv $i.sort $i
    done
    ## sort the SNP file
    for i in *.SNV.*.xls
    do
        perl $script_path/sort.variantReport.pl -i $i -o $i.sort -f Position
        mv $i.sort $i
    done
    echo `date`
fi	
