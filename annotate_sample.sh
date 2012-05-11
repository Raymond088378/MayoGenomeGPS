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
	variant_type=$( cat $run_info | grep -w '^VARIANT_TYPE' | cut -d '=' -f2)
    ## jar script to add IGV, pathway, TIssue specificity and Gene Card link
    $java/java -Xmx2g -Xms512m -jar $script_path/exome_annot.jar annotate $tool_info $output_dir/Reports_per_Sample/ $genome_version
	
    if [ $variant_type == "BOTH" -o $variant_type == "INDEL" ]
	then
		## sort the INDEL file
		for i in *.INDEL.*.xls
		do
			perl $script_path/sort.variantReport.pl -i $i -o $i.sort -f Start
			mv $i.sort $i
		done
    fi
	if [ $variant_type == "BOTH" -o $variant_type == "SNV" ]
	then
		## sort the SNP file
		for i in *.SNV.*.xls
		do
			perl $script_path/sort.variantReport.pl -i $i -o $i.sort -f Position
			mv $i.sort $i
		done
    fi
	### merge per sample files to make merged report to be uploaded to TBB
	##Merge the unfiltered file
	cd $output_dir/Reports_per_Sample/
	## SNV
	ls *.SNV.cleaned_annot.xls | sort > list
	perl $script_path/union.snv.pl list $output_dir/Reports/SNV.cleaned_annot.xls
	perl $script_path/sort.variantReport.pl -i $output_dir/Reports/SNV.cleaned_annot.xls -o $output_dir/Reports/SNV.cleaned_annot.xls.sort -f Position
	mv $output_dir/Reports/SNV.cleaned_annot.xls.sort $output_dir/Reports/SNV.cleaned_annot.xls
	
	ls *.SNV.cleaned_annot_filtered.xls | sort > list
	perl $script_path/union.snv.pl list $output_dir/Reports/SNV.cleaned_annot_filtered.xls
	perl $script_path/sort.variantReport.pl -i $output_dir/Reports/SNV.cleaned_annot_filtered.xls -o $output_dir/Reports/SNV.cleaned_annot_filtered.xls.sort -f Position
	mv $output_dir/Reports/SNV.cleaned_annot_filtered.xls.sort $output_dir/Reports/SNV.cleaned_annot_filtered.xls	
	
	
	## INDEL
	ls *.INDEL.cleaned_annot.xls | sort > list
	perl $script_path/union.indel.pl list $output_dir/Reports/INDEL.cleaned_annot.xls
	perl $script_path/sort.variantReport.pl -i $output_dir/Reports/INDEL.cleaned_annot.xls -o $output_dir/Reports/INDEL.cleaned_annot.xls.sort -f Start
	mv $output_dir/Reports/INDEL.cleaned_annot.xls.sort $output_dir/Reports/INDEL.cleaned_annot.xls
	
	ls *.INDEL.cleaned_annot_filtered.xls | sort > list
	perl $script_path/union.indel.pl list $output_dir/Reports/INDEL.cleaned_annot_filtered.xls
	perl $script_path/sort.variantReport.pl -i $output_dir/Reports/INDEL.cleaned_annot_filtered.xls -o $output_dir/Reports/INDEL.cleaned_annot_filtered.xls.sort -f Start
	mv $output_dir/Reports/INDEL.cleaned_annot_filtered.xls.sort $output_dir/Reports/INDEL.cleaned_annot_filtered.xls	
	rm list
	echo `date`
fi	
