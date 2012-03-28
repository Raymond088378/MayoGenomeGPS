#!/bin/sh

if [ $# != 4 ];
then
	echo -e "Usage: script to merge the per chr report \n<output_dir> <TempReports> <run_info> ";
else
	set -x
	echo `date`
	output_dir=$1
	TempReports=$2
	sample=$3
	run_info=$4
	tool_info=$( cat $run_info | grep -w '^TOOL_INFO' | cut -d '=' -f2)
	analysis=$( cat $run_info | grep -w '^ANALYSIS' | cut -d '=' -f2| tr "[A-Z]" "[a-z]")
    run_num=$( cat $run_info | grep -w '^OUTPUT_FOLDER' | cut -d '=' -f2)
    flowcell=`echo $run_num | awk -F'_' '{print $NF}' | sed 's/.\(.*\)/\1/'`
	java=$( cat $tool_info | grep -w '^JAVA' | cut -d '=' -f2)
    version=$( cat $run_info | grep -w '^VERSION' | cut -d '=' -f2)
	script_path=$( cat $tool_info | grep -w '^WHOLEGENOME_PATH' | cut -d '=' -f2 )

	chrs=$( cat $run_info | grep -w '^CHRINDEX' | cut -d '=' -f2)
	chrIndexes=$( echo $chrs | tr ":" "\n" )
	variant_type=$( cat $run_info | grep -w '^VARIANT_TYPE' | cut -d '=' -f2| tr "[a-z]" "[A-Z]")
	i=1
	for chr in $chrIndexes
	do
		chrArray[$i]=$chr
		let i=i+1
	done
	
	
	if [ $variant_type == "BOTH" -o $variant_type == "SNV" ]
	then
		touch $output_dir/Reports_per_Sample/$sample.SNV.report
		cat $TempReports/$sample.chr${chrArray[1]}.SNV.report >> $output_dir/Reports_per_Sample/$sample.SNV.report
		
		if [ ${#chrArray[@]} -gt 1 ]
		then
			for j in $(seq 2 ${#chrArray[@]})
			do
				cat $TempReports/$sample.chr${chrArray[$j]}.SNV.report | awk 'NR>2' >> $output_dir/Reports_per_Sample/$sample.SNV.report
			done
		fi
	fi	
	
	if [ $variant_type == "BOTH" -o $variant_type == "INDEL" ]
	then
		touch $output_dir/Reports_per_Sample/$sample.INDEL.report
		cat $TempReports/$sample.chr${chrArray[1]}.INDEL.report >> $output_dir/Reports_per_Sample/$sample.INDEL.report

		if [ ${#chrArray[@]} -gt 1 ]
		then
			for j in $(seq 2 ${#chrArray[@]})
			do
				cat $TempReports/$sample.chr${chrArray[$j]}.INDEL.report | awk 'NR>2' >> $output_dir/Reports_per_Sample/$sample.INDEL.report
			done
		fi
	fi
	
	### update the dash board
	
	if [ $analysis == "mayo" -o $analysis == "realign-mayo" ]
	then
		pos=$( cat $run_info | grep -w '^SAMPLENAMES' | cut -d '=' -f2 | tr ":" "\n" | grep -n $sample | cut -d ":" -f1)
		lanes=$( cat $run_info | grep -w '^LANEINDEX' | cut -d '=' -f2 | tr ":" "\n" | head -n $pos | tail -n 1 | tr "," " ")
		i=1
		for lane in $lanes
		do
			index=$( cat $run_info | grep -w '^LABINDEXES' | cut -d '=' -f2 | tr ":" "\n" | head -n $pos | tail -n 1 | tr "," "\n" | head -n $i | tail -n 1)
			if [ $index == "-" ]
			then
				$java/java -jar $script_path/AddSecondaryAnalysis.jar -p $script_path/AddSecondaryAnalysis.properties -l $lane -c -f $flowcell -r $run_num -s Annotation -a WholeGenome -v $version
			else
				$java/java -jar $script_path/AddSecondaryAnalysis.jar -p $script_path/AddSecondaryAnalysis.properties -l $lane -c -f $flowcell -i $index -r $run_num -s Annotation -a WholeGenome -v $version
			fi
			let i=i+1
		done		
	fi   
	echo `date`
fi	

