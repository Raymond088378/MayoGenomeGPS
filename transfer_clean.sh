#!/bin/sh

if [ $# != 2 ]
then
    echo -e "Usage: wrapper to clean intermediate files and tansfer the data to tertiary, delivery folder \n <secondary folder> < run_info >"
else
    set -x
    echo `date`
    secondary=$1
    run_info=$2
    delivery=$( cat $run_info | grep -w '^DELIVERY_FOLDER' | cut -d '=' -f2)
    tertiary=$( cat $run_info | grep -w '^TERTIARY_FOLDER' | cut -d '=' -f2)
    type=$( cat $run_info | grep -w '^TYPE' | cut -d '=' -f2 |tr "[A-Z]" "[a-z]")
    tool_info=$( cat $run_info | grep -w '^TOOL_INFO' | cut -d '=' -f2)
    sites=$( cat $tool_info | grep -w '^EMIT_ALL_SITES' | cut -d '=' -f2 | tr "[A-Z]" "[a-z]")
	multi=$( cat $run_info | grep -w '^MULTISAMPLE' | cut -d '=' -f2)
    if [ ! -s $run_info ]
    then
		echo "Runinfo file doesn't exist"
		exit 1;
    fi
    
    if [ $tertiary == "NA" ]
    then
		echo "Runinfo file doesn't have tertiary path defined"
		exit 1;
    fi
    
    if [ $delivery == "NA" ]
    then
		echo "Runinfo file doesn't have delivery path defined"
		exit 1;
    fi	
    
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
    chmod -Rf 777 $delivery/
	mv $secondary/*.html $delivery/
    if [ ! -s $delivery/Main_Document.html ]
    then
		echo "User doesn't have access to the $delivery delivery folder "
		exit 1;
    fi 
    
    if [ -d $secondary/Reports_per_Sample/temp ]
	then
		rm -R $secondary/Reports_per_Sample/temp
    fi
	
	if [ -d $secondary/Reports_per_Sample/plot ]
	then
		rm -R $secondary/Reports_per_Sample/plot
    fi
	
	mkdir $delivery/Reports_per_Sample
    chmod -Rf 777 $delivery/Reports_per_Sample
	mkdir $delivery/Reports_per_Sample/ANNOT
    chmod -Rf 777 $delivery/Reports_per_Sample/ANNOT
	if [ $type == "whole_genome" ]
    then
		mkdir $delivery/circos
		chmod -Rf 777 $delivery/circos
		mv $secondary/circos/* $delivery/circos
		mkdir $delivery/Reports_per_Sample/SV
		chmod -Rf 777 $delivery/Reports_per_Sample/SV
		mv $secondary/Reports_per_Sample/SV/*.vcf $delivery/Reports_per_Sample/SV/
		mv $secondary/Reports_per_Sample/SV/*.vcf.idx $delivery/Reports_per_Sample/SV/
		rm -R $secondary/struct
		rm -R $secondary/cnv
    fi
    
    
    mv $secondary/Reports_per_Sample/*.xls $delivery/Reports_per_Sample/
    mv $secondary/Reports_per_Sample/ANNOT/*.txt $delivery/Reports_per_Sample/ANNOT/
    mv $secondary/Reports_per_Sample/*.filter.vcf $delivery/Reports_per_Sample/
    mv $secondary/Reports_per_Sample/*.multi.vcf $delivery/Reports_per_Sample/
    
    
    rm -R $secondary/Reports_per_Sample/
    
    mkdir $delivery/Reports/
    chmod -Rf 777 $delivery/Reports/
	mv $secondary/Reports/*.xls $delivery/Reports/
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
	
	if [ $multi == "YES" ]
	then
		mv $secondary/SampleStatistics.pair.tsv $delivery/
	fi	
	
    mv $secondary/ColumnDescription_Reports.xls $delivery/
    ### make tar balls
    cd $secondary
    tar -cvzf logs.tar.gz logs
    rm -R $secondary/logs
    tar -cvzf numbers.tar.gz numbers
    rm -R $secondary/numbers
    
    ##### transfer files to tertiary folder
    mkdir -p $tertiary/variants
    if [ $sites == "yes" ]
    then
        mv $secondary/variants/*.gz $tertiary/variants/
        mv $secondary/variants/*.tbi $tertiary/variants/
    fi
    
    if [ -d $tertiary/variants ]
    then
        rm -R $secondary/variants
    fi	
    
    ### delete intermediate files
    rm -R $secondary/alignment
    rm -R $secondary/annotation
    rm -R $secondary/OnTarget
    rm -R $secondary/realign
    rm -R $secondary/TempReports
    
    if [ $delivery != "NA" ]
    then
        rm -R $secondary/IGV_BAM
    fi
    
    echo "data is transfered and intermediate files are deleted"
    echo "User needs to transfer the data to the windows share"
    echo `date`
fi	








	





