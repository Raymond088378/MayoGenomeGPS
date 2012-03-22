#!/bin/sh

if [ $# != 3 ];
then
	echo "usage:<TempReports><output_dir><run_info>";
else
	set -x
	echo `date`
	TempReports=$1
	output_dir=$2
	run_info=$3
	
	tool_info=$( cat $run_info | grep -w '^TOOL_INFO' | cut -d '=' -f2)	
	ref_flat=$( cat $tool_info | grep -w '^UCSC_REF_FLAT' | cut -d '=' -f2)
	java=$( cat $tool_info | grep -w '^JAVA' | cut -d '=' -f2 )
	variant_type=$( cat $run_info | grep -w '^VARIANT_TYPE' | cut -d '=' -f2)
	script_path=$( cat $tool_info | grep -w '^WHOLEGENOME_PATH' | cut -d '=' -f2)
	chrs=$( cat $run_info | grep -w '^CHRINDEX' | cut -d '=' -f2)
	chrIndexes=$( echo $chrs | tr ":" "\n" )
	i=1
	for chr in $chrIndexes
	do
		chrArray[$i]=$chr
		let i=i+1
	done
	
	# ## to get variant distance for snvs and indels
	for i in $(seq 1 ${#chrArray[@]})
	do 
		if [ $variant_type == "BOTH" -o $variant_type == "SNV" ]
		then
			cat $TempReports/list.chr${chrArray[$i]}.snvs | cut -f 1,2 > $TempReports/snv.chr${chrArray[$i]}.variantLocations.txt
			$java/java -Xmx2g -Xms512m -jar $script_path/exonvariantlocation.jar $ref_flat $TempReports/snv.chr${chrArray[$i]}.variantLocations.txt snp
		fi
		if [ $variant_type == "BOTH" -o $variant_type == "INDEL" ]
		then
			cat $TempReports/list.chr${chrArray[$i]}.indels | cut -f 1,2,3 > $TempReports/indel.chr${chrArray[$i]}.variantLocations.txt
			$java/java -Xmx2g -Xms512m -jar $script_path/exonvariantlocation.jar $ref_flat $TempReports/indel.chr${chrArray[$i]}.variantLocations.txt indel
		fi
	done
	if [ $variant_type == "BOTH" -o $variant_type == "SNV" ]
	then
		touch $output_dir/Reports/variantLocation_SNVs
		cat $TempReports/snv.chr${chrArray[1]}.variantLocations_out.txt >> $output_dir/Reports/variantLocation_SNVs
	fi
	if [ $variant_type == "BOTH" -o $variant_type == "IDNEL" ]
	then
		touch $output_dir/Reports/variantLocation_INDELs
		cat $TempReports/indel.chr${chrArray[1]}.variantLocations_out.txt >> $output_dir/Reports/variantLocation_INDELs	
	fi
	
	for i in $(seq 2 ${#chrArray[@]})
	do
		if [ $variant_type == "BOTH" -o $variant_type == "SNV" ]
		then
			cat $TempReports/snv.chr${chrArray[$i]}.variantLocations_out.txt | awk 'NR>8' >> $output_dir/Reports/variantLocation_SNVs
		fi
		if [ $variant_type == "BOTH" -o $variant_type == "INDEL" ]
		then
			cat $TempReports/indel.chr${chrArray[$i]}.variantLocations_out.txt | awk 'NR>8' >> $output_dir/Reports/variantLocation_INDELs
		fi
	done	
	echo `date`
fi	