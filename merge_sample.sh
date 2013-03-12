#!/bin/bash

if [ $# != 2 ]
then
    echo -e "script to merge the per sample report\nUsage : <output_dir> <run_info>";
else
	set -x
    echo `date`
    output_dir=$1
	
    run_info=$2
    tool_info=$( cat $run_info | grep -w '^TOOL_INFO' | cut -d '=' -f2)
    script_path=$( cat $tool_info | grep -w '^WORKFLOW_PATH' | cut -d '=' -f2 )
    java=$( cat $tool_info | grep -w '^JAVA' | cut -d '=' -f2 )
    genome_version=$(cat $run_info | grep -w '^GENOMEBUILD' | cut -d '=' -f2)
    variant_type=$( cat $run_info | grep -w '^VARIANT_TYPE' | cut -d '=' -f2)
    samples=$( cat $run_info | grep -w '^SAMPLENAMES' | cut -d '=' -f2 | tr ":" " ")
    multi_sample=$( cat $run_info | grep -w '^MULTISAMPLE' | cut -d '=' -f2| tr "[a-z]" "[A-Z]")
    sample_info=$( cat $run_info | grep -w '^SAMPLE_INFO' | cut -d '=' -f2)
    groups=$( cat $run_info | grep -w '^GROUPNAMES' | cut -d '=' -f2 | tr ":" " ")
    variant_type=$( cat $run_info | grep -w '^VARIANT_TYPE' | cut -d '=' -f2 |tr "[a-z]" "[A-Z]")
    snv_caller=$( cat $run_info | grep -w '^SNV_CALLER' | cut -d '=' -f2)
    ### Beauty Specific
    BEAUTYDIR="$( cat $tool_info | grep -w '^WORKFLOW_PATH' | cut -d '=' -f2 )/beauty_annot_module"
    BEAUTYDB=$( cat $tool_info | grep -w '^ANNOTATION_MODULE_DATA' | cut -d '=' -f2 )
	somatic_calling=$( cat $tool_info | grep -w '^SOMATIC_CALLING' | cut -d '=' -f2 )
	chrs=$( cat $run_info | grep -w '^CHRINDEX' | cut -d '=' -f2)
    ### merge per sample files to make merged report to be uploaded to TBB
    ##Merge the unfiltered file
    cd $output_dir/TempReports/
    mkdir -p $output_dir/Reports/
    ## SNV
    if [ $multi_sample == "NO" ]
	then
		for chr in `echo $chrs | tr ":" "\n"`
		do
			for sample in $samples
			do
				if [ $variant_type == "BOTH" -o $variant_type == "SNV" ]
				then
					if [ ! -s $sample.chr$chr.SNV.xls ]
					then
						$script_path/email.sh $sample.chr$chr.SNV.xls "doesn't exist" "sample_report.sh" $run_info
						touch $sample.chr$chr.SNV.xls.fix.log
						$script_path/wait.sh $sample.chr$chr.SNV.xls.fix.log 
					fi
					ls $sample.chr$chr.SNV.xls >> list.chr$chr.snv
					if [ ! -s $sample.chr$chr.final.SNV.xls ]
					then
						$script_path/email.sh $sample.chr$chr.final.SNV.xls "doesn't exist" "sample_report.sh" $run_info
						touch $sample.chr$chr.final.SNV.xls.fix.log
						$script_path/wait.sh $sample.chr$chr.final.SNV.xls.fix.log 
					fi
					ls $sample.chr$chr.final.SNV.xls >> list.chr$chr.final.snv
				fi
				if [ $variant_type == "BOTH" -o $variant_type == "INDEL" ]
				then
					if [ ! -s $sample.chr$chr.INDEL.xls ]
					then
						$script_path/email.sh $sample.chr$chr.INDEL.xls "doesn't exist" "sample_report.sh" $run_info
						touch $sample.chr$chr.INDEL.xls.fix.log
						$script_path/wait.sh $sample.chr$chr.INDEL.xls.fix.log 
					fi
					ls $sample.chr$chr.INDEL.xls >> list.chr$chr.indel
					if [ ! -s $sample.chr$chr.final.INDEL.xls ]
					then
						$script_path/email.sh $sample.chr$chr.final.INDEL.xls "doesn't exist" "sample_report.sh" $run_info
						touch $sample.chr$chr.final.INDEL.xls.fix.log
						$script_path/wait.sh $sample.chr$chr.final.INDEL.xls.fix.log 
					fi
					ls $sample.chr$chr.final.INDEL.xls >> list.chr$chr.final.indel
				fi
			done
		done	
		if [ $variant_type == "BOTH" -o $variant_type == "SNV" ]
		then
			for chr in `echo $chrs | tr ":" "\n" `
			do
				perl $script_path/union.snv.pl list.chr$chr.snv single $output_dir/Reports/SNV.chr$chr.xls
				perl $script_path/union.snv.pl list.chr$chr.final.snv single $output_dir/Reports/SNV.chr$chr.final.xls
				rm list.chr$chr.snv list.chr$chr.final.snv
				###raw file
				cat $output_dir/Reports/SNV.chr$chr.xls | head -2 > $output_dir/Reports/SNV.xls.header
				cat $output_dir/Reports/SNV.chr$chr.xls  | awk 'NR>2' >> $output_dir/Reports/SNV.xls
				rm $output_dir/Reports/SNV.chr$chr.xls
				### final file
				cat $output_dir/Reports/SNV.chr$chr.final.xls | head -2 > $output_dir/Reports/SNV.final.xls.header
				cat $output_dir/Reports/SNV.chr$chr.final.xls  | awk 'NR>2' >> $output_dir/Reports/SNV.final.xls
				rm $output_dir/Reports/SNV.chr$chr.final.xls
			done
			###raw file
			cat $output_dir/Reports/SNV.xls.header $output_dir/Reports/SNV.xls > $output_dir/Reports/SNV.xls.tmp
			mv $output_dir/Reports/SNV.xls.tmp $output_dir/Reports/SNV.xls
			rm $output_dir/Reports/SNV.xls.header
			### final file
			cat $output_dir/Reports/SNV.final.xls.header $output_dir/Reports/SNV.final.xls > $output_dir/Reports/SNV.final.xls.tmp
			mv $output_dir/Reports/SNV.final.xls.tmp $output_dir/Reports/SNV.final.xls
			rm $output_dir/Reports/SNV.final.xls.header
		fi
		if [ $variant_type == "BOTH" -o $variant_type == "INDEL" ]
		then
			for chr in `echo $chrs | tr ":" "\n" `
			do
				perl $script_path/union.indel.pl list.chr$chr.indel single $output_dir/Reports/INDEL.chr$chr.xls
				perl $script_path/union.indel.pl list.final.chr$chr.indel single $output_dir/Reports/INDEL.chr$chr.final.xls
				rm list.chr$chr.indel list.chr$chr.final.indel
				###raw file
				cat $output_dir/Reports/INDEL.chr$chr.xls | head -2 > $output_dir/Reports/INDEL.xls.header
				cat $output_dir/Reports/INDEL.chr$chr.xls  | awk 'NR>2' >> $output_dir/Reports/INDEL.xls
				rm $output_dir/Reports/INDEL.chr$chr.xls
				### final file
				cat $output_dir/Reports/INDEL.chr$chr.final.xls | head -2 > $output_dir/Reports/INDEL.final.xls.header
				cat $output_dir/Reports/INDEL.chr$chr.final.xls  | awk 'NR>2' >> $output_dir/Reports/INDEL.final.xls
				rm $output_dir/Reports/INDEL.chr$chr.final.xls
			done
			###raw file
			cat $output_dir/Reports/INDEL.xls.header $output_dir/Reports/INDEL.xls > $output_dir/Reports/INDEL.xls.tmp
			mv $output_dir/Reports/INDEL.xls.tmp $output_dir/Reports/INDEL.xls
			rm $output_dir/Reports/INDEL.xls.header
			### final file
			cat $output_dir/Reports/INDEL.final.xls.header $output_dir/Reports/INDEL.final.xls > $output_dir/Reports/INDEL.final.xls.tmp
			mv $output_dir/Reports/INDEL.final.xls.tmp $output_dir/Reports/INDEL.final.xls
			rm $output_dir/Reports/INDEL.final.xls.header
		fi
	else
		for chr in `echo $chrs | tr ":" "\n" `
		do
			for group in $groups
			do
				if [ ! -s $group.chr$chr.SNV.xls ]
				then
					$script_path/email.sh $group.chr$chr.SNV.xls "doesn't exist" "sample_report.sh" $run_info
					touch $group.chr$chr.SNV.xls.fix.log
					$script_path/wait.sh $group.chr$chr.SNV.xls.fix.log 
				fi
				ls $group.chr$chr.SNV.xls >> list.chr$chr.snv
				if [ ! -s $group.chr$chr.final.SNV.xls ]
				then
					$script_path/email.sh $group.chr$chr.final.SNV.xls "doesn't exist" "sample_report.sh" $run_info
					touch $group.chr$chr.final.SNV.xls.fix.log
					$script_path/wait.sh $group.chr$chr.final.SNV.xls.fix.log 
				fi
				ls $group.chr$chr.final.SNV.xls >> list.chr$chr.final.snv
				if [ ! -s $group.chr$chr.INDEL.xls ]
				then
					$script_path/email.sh $group.chr$chr.INDEL.xls "doesn't exist" "sample_report.sh" $run_info
					touch $group.chr$chr.INDEL.xls.fix.log
					$script_path/wait.sh $group.chr$chr.INDEL.xls.fix.log 
				fi
				ls $group.chr$chr.INDEL.xls >> list.chr$chr.indel
				if [ ! -s $group.chr$chr.final.INDEL.xls ]
				then
					$script_path/email.sh $group.chr$chr.final.INDEL.xls "doesn't exist" "sample_report.sh" $run_info
					touch $group.chr$chr.final.INDEL.xls.fix.log
					$script_path/wait.sh $group.chr$chr.final.INDEL.xls.fix.log 
				fi
				ls $group.chr$chr.final.INDEL.xls >> list.chr$chr.final.indel
			done
		done
		for chr in `echo $chrs | tr ":" "\n"`
		do
			perl $script_path/union.snv.pl list.chr$chr.snv multi $output_dir/Reports/SNV.chr$chr.xls
			perl $script_path/union.snv.pl list.chr$chr.final.snv multi $output_dir/Reports/SNV.chr$chr.final.xls
			perl $script_path/union.indel.pl list.chr$chr.indel multi $output_dir/Reports/INDEL.chr$chr.xls
			perl $script_path/union.indel.pl list.chr$chr.final.indel multi $output_dir/Reports/INDEL.chr$chr.final.xls	
			rm list.chr$chr.snv list.chr$chr.final.snv list.chr$chr.indel list.chr$chr.final.indel	
			###raw file
			cat $output_dir/Reports/INDEL.chr$chr.xls | head -2 > $output_dir/Reports/INDEL.xls.header
			cat $output_dir/Reports/INDEL.chr$chr.xls  | awk 'NR>2' >> $output_dir/Reports/INDEL.xls
			rm $output_dir/Reports/INDEL.chr$chr.xls
			### final file
			cat $output_dir/Reports/INDEL.chr$chr.final.xls | head -2 > $output_dir/Reports/INDEL.final.xls.header
			cat $output_dir/Reports/INDEL.chr$chr.final.xls  | awk 'NR>2' >> $output_dir/Reports/INDEL.final.xls
			rm $output_dir/Reports/INDEL.chr$chr.final.xls
			###raw file
			cat $output_dir/Reports/SNV.chr$chr.xls | head -2 > $output_dir/Reports/SNV.xls.header
			cat $output_dir/Reports/SNV.chr$chr.xls  | awk 'NR>2' >> $output_dir/Reports/SNV.xls
			rm $output_dir/Reports/SNV.chr$chr.xls
			### final file
			cat $output_dir/Reports/SNV.chr$chr.final.xls | head -2 > $output_dir/Reports/SNV.final.xls.header
			cat $output_dir/Reports/SNV.chr$chr.final.xls  | awk 'NR>2' >> $output_dir/Reports/SNV.final.xls
			rm $output_dir/Reports/SNV.chr$chr.final.xls
		done
		###raw file
		cat $output_dir/Reports/INDEL.xls.header $output_dir/Reports/INDEL.xls > $output_dir/Reports/INDEL.xls.tmp
		mv $output_dir/Reports/INDEL.xls.tmp $output_dir/Reports/INDEL.xls
		rm $output_dir/Reports/INDEL.xls.header
		### final file
		cat $output_dir/Reports/INDEL.final.xls.header $output_dir/Reports/INDEL.final.xls > $output_dir/Reports/INDEL.final.xls.tmp
		mv $output_dir/Reports/INDEL.final.xls.tmp $output_dir/Reports/INDEL.final.xls
		rm $output_dir/Reports/INDEL.final.xls.header
		###raw file
		cat $output_dir/Reports/SNV.xls.header $output_dir/Reports/SNV.xls > $output_dir/Reports/SNV.xls.tmp
		mv $output_dir/Reports/SNV.xls.tmp $output_dir/Reports/SNV.xls
		rm $output_dir/Reports/SNV.xls.header
		### final file
		cat $output_dir/Reports/SNV.final.xls.header $output_dir/Reports/SNV.final.xls > $output_dir/Reports/SNV.final.xls.tmp
		mv $output_dir/Reports/SNV.final.xls.tmp $output_dir/Reports/SNV.final.xls
		rm $output_dir/Reports/SNV.final.xls.header
		### Merge the TUMOR files
		cd $output_dir/Reports_per_Sample
		if [ $somatic_calling == "YES" ]
		then
			for group in $groups
			do
				if [ ! -s TUMOR.$group.SNV.xls ]
				then
					$script_path/email.sh TUMOR.$group.SNV.xls "doesn't exist" "sample_report.sh" $run_info
					touch TUMOR.$group.SNV.xls.fix.log
					$script_path/wait.sh TUMOR.$group.SNV.xls.fix.log 
				fi
				ls TUMOR.$group.SNV.xls >> list.snv
				if [ ! -s TUMOR.$group.SNV.final.xls ]
				then
					$script_path/email.sh TUMOR.$group.SNV.final.xls "doesn't exist" "sample_report.sh" $run_info
					touch TUMOR.$group.SNV.final.xls.fix.log
					$script_path/wait.sh TUMOR.$group.SNV.final.xls.fix.log 
				fi
				ls TUMOR.$group.SNV.final.xls >> list.final.snv
				if [ ! -s TUMOR.$group.INDEL.xls ]
				then
					$script_path/email.sh TUMOR.$group.INDEL.xls "doesn't exist" "sample_report.sh" $run_info
					touch TUMOR.$group.INDEL.xls.fix.log
					$script_path/wait.sh TUMOR.$group.INDEL.xls.fix.log 
				fi
				ls TUMOR.$group.INDEL.xls >> list.indel
				if [ ! -s $group.INDEL.final.xls ]
				then
					$script_path/email.sh TUMOR.$group.INDEL.final.xls "doesn't exist" "sample_report.sh" $run_info
					touch TUMOR.$group.INDEL.final.xls.fix.log
					$script_path/wait.sh TUMOR.$group.INDEL.final.xls.fix.log 
				fi
				ls TUMOR.$group.INDEL.final.xls >> list.final.indel
			done
			perl $script_path/union.snv.pl list.snv multi $output_dir/Reports/TUMOR.SNV.xls
			perl $script_path/union.snv.pl list.final.snv multi $output_dir/Reports/TUMOR.SNV.final.xls
			perl $script_path/union.indel.pl list.indel multi $output_dir/Reports/TUMOR.INDEL.xls
			perl $script_path/union.indel.pl list.final.indel multi $output_dir/Reports/TUMOR.INDEL.final.xls	
			rm list.snv list.final.snv list.indel list.final.indel
		fi	
	fi
	if [ $snv_caller == "BEAUTY_EXOME" ]
	then
	
		###Intended to reduce repitious coding.
		callJasonsScripts(){
			INFILTER="$output_dir/Reports/$1.final.xls"
			CANCTMP="$output_dir/Reports/$1.cancerDrug.xls"
			PANLTMP="$output_dir/Reports/$1.clinPanel.xls"
			KINOMTPM="$output_dir/Reports/$1.kinome.xls"

			perl $BEAUTYDIR/addCancerDrugs.pl $BEAUTYDB/NCCN.CancerDrugs.txt $INFILTER $CANCTMP
			perl $BEAUTYDIR/addClinPanel.pl $BEAUTYDB/uniq_genes_in_clinical_panels.txt $CANCTMP $PANLTMP
			perl $BEAUTYDIR/addKinome.pl $BEAUTYDB/kinome.genes.txt $PANLTMP $KINOMTPM

			rm $INFILTER $CANCTMP $PANLTMP
			mv $KINOMTPM $INFILTER
		}
	fi

	### Adding Beauty Anntation Module
	if [ $snv_caller == "BEAUTY_EXOME" ]
	then
		callJasonsScripts "SNV"
		callJasonsScripts "INDEL"
		## Script to create additional annotation file, with pharmicogentic emphasis
		$script_path/PharmacoAnnotModule.sh $output_dir/Reports/SNV.final.xls $run_info
		$script_path/PharmacoAnnotModule.sh $output_dir/Reports/INDEL.final.xls $run_info
		
		if [ $multi_sample != "NO" ] ; then
			callJasonsScripts "TUMOR.SNV"
			callJasonsScripts "TUMOR.INDEL"
			$script_path/PharmacoAnnotModule.sh $output_dir/Reports/TUMOR.SNV.final.xls $run_info
			$script_path/PharmacoAnnotModule.sh $output_dir/Reports/TUMOR.final.INDEL.xls $run_info
		fi
	fi

    echo `date`
fi	
