#!/bin/sh

if [ $# != 3 ];
then
    echo "Usage: <output_dir><TempReports><run_info>";
else
    set -x
    echo `date`
    output_dir=$1
    TempReports=$2
    run_info=$3
    
    tool_info=$( cat $run_info | grep -w '^TOOL_INFO' | cut -d '=' -f2)	
    script_path=$( cat $tool_info | grep -w '^WHOLEGENOME_PATH' | cut -d '=' -f2 )
    variant_type=$( cat $run_info | grep -w '^VARIANT_TYPE' | cut -d '=' -f2| tr "[a-z]" "[A-Z]")
	java=$( cat $tool_info | grep -w '^JAVA' | cut -d '=' -f2 )
    genome_version=$(cat $run_info | grep -w '^GENOMEBUILD' | cut -d '=' -f2)
	chrs=$( cat $run_info | grep -w '^CHRINDEX' | cut -d '=' -f2)
    chrIndexes=$( echo $chrs | tr ":" "\n" )
    i=1
    for chr in $chrIndexes
    do
            chrArray[$i]=$chr
            let i=i+1
    done
    
    if [ $variant_type == "BOTH" -o $variant_type == "SNV" ]
	then
		touch $output_dir/Reports/SNV.report
		cat $TempReports/list.chr${chrArray[1]}.snvs.SNV.report >> $output_dir/Reports/SNV.report
	fi
	
	if [ $variant_type == "BOTH" -o $variant_type == "INDEL" ]
	then
		touch $output_dir/Reports/INDEL.report
		cat $TempReports/list.chr${chrArray[1]}.indels.INDEL.report >> $output_dir/Reports/INDEL.report
    fi
	
    if [ ${#chrArray[@]} -gt 1 ]
	then
		for j in $(seq 2 ${#chrArray[@]})
		do
			if [ $variant_type == "BOTH" -o $variant_type == "SNV" ]
			then
				cat $TempReports/list.chr${chrArray[$j]}.snvs.SNV.report | awk 'NR>2' >> $output_dir/Reports/SNV.report
			fi
			if [ $variant_type == "BOTH" -o $variant_type == "INDEL" ]
			then	
				cat $TempReports/list.chr${chrArray[$j]}.indels.INDEL.report | awk 'NR>2' >> $output_dir/Reports/INDEL.report
			fi	
		done
    fi
    $java/java -Xmx2g -Xms512m -jar $script_path/exome_annot.jar annotate $tool_info $output_dir/Reports/ $genome_version
    cd $output_dir/Reports/
    
	## sort the indel files
	if [ $variant_type == "BOTH" -o $variant_type == "INDEL" ]
	then
		for i in INDEL.*.xls
		do
			perl $script_path/sort.variantReport.pl -i $i -o $i.sort -f Start
			mv $i.sort $i
		done
    fi
	
	## sort the SNP file
    if [ $variant_type == "BOTH" -o $variant_type == "SNV" ]
	then
		for i in SNV.*.xls
		do
			perl $script_path/sort.variantReport.pl -i $i -o $i.sort -f Position
			mv $i.sort $i
		done
    fi
	echo `date`
fi	


