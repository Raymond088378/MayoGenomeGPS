#!/bin/bash
#	INFO	
#	script generates HTML report and coverage plot graph

if [ $# != 2 ]
then
	echo -e "script to generate html and send an email to the user of workflow completion\nUsage: <output dir> <run info file>";
else
	set -x
	echo `date`
	output_dir=$1
	run_info=$2
	
	tool_info=$( cat $run_info | grep -w '^TOOL_INFO' | cut -d '=' -f2)
	memory_info=$( cat $run_info | grep -w '^MEMORY_INFO' | cut -d '=' -f2)
	PI=$( cat $run_info | grep -w '^PI' | cut -d '=' -f2)
	run_num=$( cat $run_info | grep -w '^OUTPUT_FOLDER' | cut -d '=' -f2)
	analysis=$( cat $run_info | grep -w '^ANALYSIS' | cut -d '=' -f2| tr "[A-Z]" "[a-z]")
	script_path=$( cat $tool_info | grep -w '^WORKFLOW_PATH' | cut -d '=' -f2)
	java=$( cat $tool_info | grep -w '^JAVA' | cut -d '=' -f2 )
    run_num=$( cat $run_info | grep -w '^OUTPUT_FOLDER' | cut -d '=' -f2)
	type=$( cat $run_info | grep -w '^TYPE' | cut -d '=' -f2|tr "[a-z]" "[A-Z]")
	upload_tb=$( cat $tool_info | grep -w '^UPLOAD_TABLEBROWSER' | cut -d '=' -f2| tr "[a-z]" "[A-Z]")
	tool=$( cat $run_info | grep -w '^TYPE' | cut -d '=' -f2| tr "[A-Z]" "[a-z]")
    samples=$( cat $run_info | grep -w '^SAMPLENAMES' | cut -d '=' -f2 | tr ":" " ")    
	flowcell=`echo $run_num | awk -F'_' '{print $NF}' | sed 's/.\(.*\)/\1/'`
	# generate Coverage plot
	for sample in $samples
	do
		$script_path/dashboard.sh $sample $run_info Results started
	done
	cd $output_dir/numbers
	
	if [[ $analysis != "alignment" && $analysis != "annotation"  && $analysis != "ontarget" ]] 
	then
		$script_path/generate.coverage.sh $output_dir/numbers $output_dir $run_info
	fi
	
	if [[ $analysis != "alignment" && $analysis != "annotation"  && $analysis != "ontarget" ]] 
	then
		$script_path/create.igv.pl -o $output_dir -r $run_info
	fi
	
	$script_path/MainDocument.pl -r $run_info -p $output_dir
	### generate readme file
	$script_path/generate_readme.sh $output_dir $run_info
	END=`date`
	echo -e "Analysis Ends at :" >> $output_dir/log.txt
	echo -e "${END}" >>  $output_dir/log.txt
	cd $output_dir/logs
	
	files=`ls -lhrt | awk 'NR>1' |awk -F' ' '{print $NF}' | tr "\n" " "`
	cat $files | grep -w '^ERROR ' > $output_dir/errorlog
	cat  $files | grep -w '^WARNING ' > $output_dir/warninglog	
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
	
	email=`finger $USER | awk -F ';' '{print $2}'`
	SUB="$tool workflow completion for RunID ${run_num} "
	MESG=" ${text} ${text1} $tool workflow completed for ${run_num} on ${END} and ready for tertiary analysis in ${output_dir} "
	## send the completion email
	TO=$USER
	echo -e "$MESG\n\nTIMESTAMPS:" | cat - $output_dir/log.txt | mailx -s "$SUB" -c Kahl.Jane@mayo.edu,Hossain.Asif@mayo.edu "$email" 
	for sample in $samples
	do
		$script_path/dashboard.sh $sample $run_info Results complete
	done
	if [ $type == "exome" ]
	then
		tool=Exome
	else
		tool=WholeGenome
    fi
	if [[ $analysis == "mayo" || $analysis == "realign-mayo" ]]
	then
		mem=$( cat $memory_info | grep -w '^AddSecondaryAnalysis_JVM' | cut -d '=' -f2)
		$java/java $mem -jar $script_path/AddSecondaryAnalysis.jar -p $script_path/AddSecondaryAnalysis.properties -c -f $flowcell -r $run_num -s Complete -a $tool
	fi
	echo `date`
fi

	
	
