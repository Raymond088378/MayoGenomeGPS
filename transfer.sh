#!/bin/sh

if [ $# != 4 ]
then
	echo -e "Usage: wrapper to clean intermediate files and tansfer the data to tertiary, delivery folder \n <secondary folder> <delivery folder><tertiary folder> <whole_genome/exome>"
else
	set -x
	echo `date`
	secondary=$1
	delivery=$2
	tertiary=$3
	type=$4
	type=`echo $type | tr "[A-Z]" "[a-z]"`
	
	if [ ! -d $secondary ]
	then
		echo " $secondary secondary folder doesn't exist"
		exit 1;
	fi
	
	if [ ! -d $delivery ]
	then
		echo " $delivery delivery folder doesn't exist"
		exit 1;
	fi

	if [ ! -d $tertiary ]
	then
		mkdir -p $tertiary
		echo "$tertiary tertiary folder created"
	fi
		
	### transfer the data to delivery folder
	mv $secondary/*.html $delivery/
	if [ ! -s $delivery/Main_Document.html ]
	then
		echo "User doesn't have access to the $delivery delivery folder "
		exit 1;
	fi 
	
	rm $secondary/Reports_per_Sample/*.report
	rm -R $secondary/Reports_per_Sample/temp
	rm -R $secondary/Reports_per_Sample/plot
	mkdir $delivery/Reports_per_Sample
	mkdir $delivery/Reports_per_Sample/ANNOT
	if [ $type == "whole_genome" ]
	then
		mkdir $delivery/Reports_per_Sample/SV
		mv $secondary/Reports_per_Sample/SV/*.vcf $delivery/Reports_per_Sample/SV/
		mv $secondary/Reports_per_Sample/SV/*.vcf.idx $delivery/Reports_per_Sample/SV/
	fi
	
	mv $secondary/Reports_per_Sample/*.xls $delivery/Reports_per_Sample
	mv $secondary/Reports_per_Sample/ANNOT/*.txt $delivery/Reports_per_Sample/ANNOT/
	mv $secondary/Reports_per_Sample/*.filter.vcf $delivery/Reports_per_Sample/
	mv $secondary/Reports_per_Sample/*.filter.vcf.idx $delivery/Reports_per_Sample/
	mkdir $tertiary/Reports_per_Sample/
	mv $secondary/Reports_per_Sample/*.raw.vcf $tertiary/Reports_per_Sample/
	mv $secondary/Reports_per_Sample/*.raw.vcf.idx $tertiary/Reports_per_Sample/	
	
	rm -R $secondary/Reports_per_Sample/
	
	rm $secondary/Reports/*.report
	mkdir $delivery/Reports/
	mv $secondary/Reports/* $delivery/Reports/
	rm -R $secondary/Reports/
	mv $secondary/Coverage.JPG $delivery/	
	if [ $type == "exome" ]
	then
		mv $secondary/exome_workflow.png $delivery/
	else
		mv $secondary/whole_genome.png $delivery/
	fi
	mv $secondary/igv_session.xml $delivery/
	mv $secondary/IGV_Setup.doc $delivery/
	mv $secondary/SampleStatistics.tsv $delivery/
	
	### make tar balls
	cd $secondary
	tar -cvzf logs.tar.gz logs
	rm -R $secondary/logs
	tar -cvzf numbers.tar.gz numbers
	rm -R $secondary/numbers
	
	##### transfer files to tertiary folder
	cp -R $secondary/variants $tertiary/
	if [ -d $tertiary/variants ]
	then
		rm -R $secondary/variants
	fi	
	
	### delete intermediate files
	rm -R $secondary/alignment
	rm -R $secondary/annotation
	rm -R $secondary/job_ids
	rm -R $secondary/OnTarget
	rm -R $secondary/realign
	rm -R $secondary/TempReports
	
	echo "data is transfered and intermediate files are deleted"	
	echo `date`
fi	








	





