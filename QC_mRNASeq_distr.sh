#!/bin/sh

if [ $# != 4 ]
then
	echo "Usage:";
	echo "1. if more than one flowcell, provide full path to text file listing flowcells --> just the flowcell names. Else provide flowcell name";
	echo "2. output directory for plot";
	echo "3. script path";
	echo "4. 0 or 1 value --> 0 if plot required for only one flowcell, 1 if plot required for more than one flowcell";
else
	set -x
	echo `date`
	list=$1
	output_dir=$2
	script_path=$3
	flag=$4
	
	if [ $flag == 1 ]
	then
		echo "More than one flowcell"
		for i in `cat $list`
		do
			current=$output_dir/$i
			for j in 1 3 4 5 6
			do 
				##extracting sample names, used reads, mapped reads, genome mapped and jn mapped
				cat $current/SampleStatistics.tsv | sed -n "$j p" | tr "\t" "\n" | sed 's/,//g' > $current/$j.$i.txt 
			done
			for j in 1 3 4 5 6
			do 
				sed -i '1d' $current/$j.$i.txt
			done
			for j in 4 5 6
			do 
				## removing % values
				cat $current/$j.$i.txt | cut -f1 -d "(" > $current/$j.$i.new.txt 
			done
			rm $current/4.$i.txt $current/5.$i.txt $current/6.$i.txt
			paste $current/3.$i.txt $current/4.$i.new.txt $current/5.$i.new.txt $current/6.$i.new.txt >> $current/plot.$i.txt
			rm $current/3.$i.txt $current/4.$i.new.txt $current/5.$i.new.txt $current/6.$i.new.txt
		done
		
		echo -e "UsedReads\tMappedReads\tGenomeMapped\tJunctionMapped" >> $output_dir/final_plot.tmp
		# echo -e "Samples" >> $output_dir/final_samples.tmp
		for i in `cat $list`
		do
			current=$output_dir/$i
			cat $current/plot.$i.txt >> $output_dir/final_plot.tmp
			cat $current/1.$i.txt >> $output_dir/final_samples.tmp
			rm $current/1.$i.txt $current/plot.$i.txt
		done
		## removing empty lines
		awk NF $output_dir/final_plot.tmp > $output_dir/final_plot.txt
		awk NF $output_dir/final_samples.tmp > $output_dir/final_samples.txt
		rm $output_dir/final_plot.tmp $output_dir/final_samples.tmp 
		## plotting distribution
		Rscript $script_path/QC_mRNASeq_Rplot.r $output_dir/final_plot.txt $output_dir/final_samples.txt
		rm $output_dir/final_plot.txt $output_dir/final_samples.txt
	else
		echo "One flowcell"
		for i in $list
		do
			current=$output_dir
			cd $current
			for j in 1 3 4 5 6
			do 
				##extracting sample names, used reads, mapped reads, genome mapped and jn mapped
				cat $current/SampleStatistics.tsv | sed -n "$j p" | tr "\t" "\n" | sed 's/,//g' > $current/$j.$i.txt 
			done
			for j in 1 3 4 5 6
			do 
				sed -i '1d' $current/$j.$i.txt
			done
			for j in 4 5 6
			do 
				cat $current/$j.$i.txt | cut -f1 -d "(" > $current/$j.$i.new.txt
			done
			rm $current/4.$i.txt $current/5.$i.txt $current/6.$i.txt
			echo -e "UsedReads\tMappedReads\tGenomeMapped\tJunctionMapped" >> $current/plot.$i.tmp
			paste $current/3.$i.txt $current/4.$i.new.txt $current/5.$i.new.txt $current/6.$i.new.txt >> $current/plot.$i.tmp
			rm $current/3.$i.txt $current/4.$i.new.txt $current/5.$i.new.txt $current/6.$i.new.txt
			## removing empty lines
			awk NF $current/plot.$i.tmp > $current/plot.$i.txt
			## plotting distribution
			Rscript $script_path/QC_mRNASeq_Rplot.r $current/plot.$i.txt $current/1.$i.txt
			rm $current/plot.$i.txt $current/1.$i.txt $current/plot.$i.tmp
		done
	fi
fi