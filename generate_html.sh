#!/bin/sh
#	INFO	
#	script generates HTML report and coverage plot graph

if [ $# != 2 ]
then
	echo "Usage: <output dir> <run info file>";
else
	set -x
	echo `date`
	output_dir=$1
	run_info=$2
	
	tool_info=$( cat $run_info | grep -w '^TOOL_INFO' | cut -d '=' -f2)
	PI=$( cat $run_info | grep -w '^PI' | cut -d '=' -f2)
	run_num=$( cat $run_info | grep -w '^OUTPUT_FOLDER' | cut -d '=' -f2)
	analysis=$( cat $run_info | grep -w '^ANALYSIS' | cut -d '=' -f2)
	script_path=$( cat $tool_info | grep -w '^WHOLEGENOME_PATH' | cut -d '=' -f2)
	java=$( cat $tool_info | grep -w '^JAVA' | cut -d '=' -f2 )
	http_server=$( cat $tool_info | grep -w '^HTTP_SERVER' | cut -d '=' -f2 )
	queue=$( cat $run_info | grep -w '^QUEUE' | cut -d '=' -f2)
    run_num=$( cat $run_info | grep -w '^OUTPUT_FOLDER' | cut -d '=' -f2)
	lanes=$( cat $run_info | grep -w '^LANEINDEX' | cut -d '=' -f2)
	tool=$( cat $run_info | grep -w '^TYPE' | cut -d '=' -f2)
	type=$( cat $run_info | grep -w '^TOOL' | cut -d '=' -f2|tr "[a-z]" "[A-Z]")
	email=$( cat $run_info | grep -w '^EMAIL' | cut -d '=' -f2)
	GenomeBuild=$( cat $run_info | grep -w '^GENOMEBUILD' | cut -d '=' -f2)
	variant_type=$( cat $run_info | grep -w '^VARIANT_TYPE' | cut -d '=' -f2)
	UCSC=$( cat $tool_info | grep -w '^UCSC_TRACKS' | cut -d '=' -f2)
	analysis=`echo "$analysis" | tr "[A-Z]" "[a-z]"`
	upload_tb=$( cat $run_info | grep -w '^UPLOAD_TABLEBROWSER' | cut -d '=' -f2)
	upload_tb=`echo "$upload_tb" | tr "[a-z]" "[A-Z]"`
	variant_type=`echo "$variant_type" | tr "[a-z]" "[A-Z]"`
	tool=$( cat $run_info | grep -w '^TYPE' | cut -d '=' -f2)
	tool=`echo "$tool" | tr "[A-Z]" "[a-z]"`
	CaptureKit=$( cat $tool_info | grep -w '^CAPTUREKIT' | cut -d '=' -f2 )
	master_gene_file=$( cat $tool_info | grep -w '^MASTER_GENE_FILE' | cut -d '=' -f2 )
	samples=$( cat $run_info | grep -w '^SAMPLENAMES' | cut -d '=' -f2 | tr ":" " " )
	bed=$( cat $tool_info | grep -w '^BEDTOOLS' | cut -d '=' -f2 )
	if [ $tool == "whole_genome" ]
	then
		kit=$output_dir/bed_file.bed
	else
		kit=$CaptureKit
	fi    
	# generate Coverage plot
	
	cd $output_dir/numbers
	if [[ $analysis != "alignment" && $analysis != "annotation" ]] 
	then
		region=`awk '{sum+=$3-$2+1; print sum}' $kit | tail -1`
		Rscript $script_path/coverage_plot.r $region $samples
		mv $output_dir/numbers/coverage.jpeg $output_dir/Coverage.JPG
	fi
	#rm $output_dir/bed_file.bed
	if [ $analysis != "annotation" -a $analysis != "alignment" ]
	then
		perl $script_path/create.igv.pl -o $output_dir -r $run_info
	fi
	perl $script_path/MainDocument.pl -r $run_info -p $output_dir
	
	## create tsv file for sample statistcs
	perl $script_path/SampleStatistics.pl -r $run_info -p $output_dir
	## TableBrowser upload
	if [ $upload_tb == "YES" ]
	then
		PI_LANID=$( echo $PI | cut -d '_' -f 3 )
		if [ $GenomeBuild == "hg19" ]
		then
			$java/java -Xmx7g -Xms512m -jar $script_path/TREATUploader.jar -n $PI_LANID -u $run_num -i $output_dir/Reports/INDEL.cleaned_annot.xls -s $output_dir/Reports/SNV.cleaned_annot.xls -r $run_num -f $script_path/format.SNV.txt
			else
			$java/java -Xmx7g -Xms512m -jar $script_path/TREATUploader.jar -n $PI_LANID -u $run_num -i $output_dir/Reports/INDEL.cleaned_annot.xls -s $output_dir/Reports/SNV.cleaned_annot.xls -r $run_num
		fi		
		echo -e "Variants uploaded to TableBrowser" >> $output_dir/log.txt
	else
		echo -e "Variants Not uploaded to TableBrowser" >> $output_dir/log.txt
	fi	
	END=`date`
	echo -e "Analysis Ends at :" >> $output_dir/log.txt
	echo -e "${END}" >>  $output_dir/log.txt
	cd $output_dir
    logs=`find -name 'logs' | sed -e '/\.\//s///g'`
    for i in $logs
    do
        cat $i/$type* >> $output_dir/LOG
    done
    cat $output_dir/LOG | grep -w '^ERROR' > $output_dir/errorlog
	cat $output_dir/LOG | grep -w '^WARNING' > $output_dir/warninglog	
    rm $output_dir/LOG
	
	e_size=`ls -l $output_dir/errorlog | awk '{ print $5 }'`
	w_size=`ls -l $output_dir/warninglog | awk '{ print $5 }'`
	if [ $e_size -le 0 ]
	then
		text="SUCCESS"
	else
		text="ERROR"
	fi	
	if [ $w_size -le 0 ]
	then
		text1="with no warnings"
	else
		text1="with warnings"
	fi		
	
	TO=`id | cut -d '(' -f2 | sed -e '/)/s///g' | cut -d ' ' -f1`
	SUB="$tool workflow completion for RunID ${run_num} "
	MESG=" ${text} ${text1} $tool workflow completed for ${run_num} on ${END} and ready for tertiary analysis in ${output_dir} "
	## send the completion email
	echo -e "$MESG" | mailx -v -s "$SUB" "$TO" 
	echo `date`
fi

	
	
