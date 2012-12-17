#!/bin/sh

if [ $# != 2 ]
then
	echo "\nUsage: <list of flowcell names> </path/to/delivery folder for PI>";
else
	set -x
	echo `date`
	list=$1
	output_dir=$2
	
	for i in `cat $list`
	do
		output=$output_dir/$i/secondary
		cd $output
		mkdir $output/QC
		reports=$output/Reports
		QC=$output/QC
		
		## extracting sample specific columns to count for zeros
		
		head -1 $reports/GeneCount.tsv	| cut -f3 > $reports/samlist
		for j in `cat $reports/samlist`; do echo -e ${j}_NoReadsMapped >> $QC/$j.noreads; done
		for i in 1; do val=`expr $i + 2`; awk '{if ($'$val' == 0) print $2}' $reports/GeneCount.tsv >> $QC/$j.noreads; done
		rm $reports/samlist
		
		head -1 $reports/GeneCount.tsv	| cut -f4 > $reports/samlist
		for j in `cat $reports/samlist`; do echo -e ${j}_NoReadsMapped >> $QC/$j.noreads; done
		for i in 2; do val=`expr $i + 2`; awk '{if ($'$val' == 0) print $2}' $reports/GeneCount.tsv >> $QC/$j.noreads; done
		rm $reports/samlist
		
		head -1 $reports/GeneCount.tsv	| cut -f5 > $reports/samlist
		for j in `cat $reports/samlist`; do echo -e ${j}_NoReadsMapped >> $QC/$j.noreads; done
		for i in 3; do val=`expr $i + 2`; awk '{if ($'$val' == 0) print $2}' $reports/GeneCount.tsv >> $QC/$j.noreads; done
		rm $reports/samlist
		
		head -1 $reports/GeneCount.tsv	| cut -f6 > $reports/samlist
		for j in `cat $reports/samlist`; do echo -e ${j}_NoReadsMapped >> $QC/$j.noreads; done
		for i in 4; do val=`expr $i + 2`; awk '{if ($'$val' == 0) print $2}' $reports/GeneCount.tsv >> $QC/$j.noreads; done
		rm $reports/samlist
		
		head -1 $reports/GeneCount.tsv	| cut -f7 > $reports/samlist
		for j in `cat $reports/samlist`; do echo -e ${j}_NoReadsMapped >> $QC/$j.noreads; done
		for i in 5; do val=`expr $i + 2`; awk '{if ($'$val' == 0) print $2}' $reports/GeneCount.tsv >> $QC/$j.noreads; done
		rm $reports/samlist
		
		head -1 $reports/GeneCount.tsv	| cut -f8 > $reports/samlist
		for j in `cat $reports/samlist`; do echo -e ${j}_NoReadsMapped >> $QC/$j.noreads; done
		for i in 6; do val=`expr $i + 2`; awk '{if ($'$val' == 0) print $2}' $reports/GeneCount.tsv >> $QC/$j.noreads; done
		rm $reports/samlist
		
		head -1 $reports/GeneCount.tsv	| cut -f9 > $reports/samlist
		for j in `cat $reports/samlist`; do echo -e ${j}_NoReadsMapped >> $QC/$j.noreads; done
		for i in 7; do val=`expr $i + 2`; awk '{if ($'$val' == 0) print $2}' $reports/GeneCount.tsv >> $QC/$j.noreads; done
		rm $reports/samlist
		
		head -1 $reports/GeneCount.tsv	| cut -f10 > $reports/samlist
		for j in `cat $reports/samlist`; do echo -e ${j}_NoReadsMapped >> $QC/$j.noreads; done
		for i in 8; do val=`expr $i + 2`; awk '{if ($'$val' == 0) print $2}' $reports/GeneCount.tsv >> $QC/$j.noreads; done
		rm $reports/samlist
	done

fi