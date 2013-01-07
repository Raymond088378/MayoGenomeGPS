#!/bin/bash
	
########################################################
###### 	GENE SUMMARY FOR TUMOR/NORMAL PAIR WHOLE GENOME ANALYSIS PIPELINE

######		Program:			gene.summary.sh
######		Date:				11/09/2011
######		Summary:			Summary of SNV, INDEL, SV and CNV per gene
######		Input 
######		$1	=	output directory
######		$2	=	SNV directory
######		$3	=	INDEL directory
######		$4	=	CNV directory
######		$5	=	SV directory
######		$6	=	/path/to/run_info.txt

########################################################

if [ $# != 4 ]
then
    echo -e "script to get the gene summary file with information about the vaianats and SVs\nUsage: </path/to/output directory> <sample name> </path/to/run_info.txt> </path/yo/Reports_per_Sample>";
else
    set -x
    echo `date`
    output_dir=$1
    sample=$2
    run_info=$3
    output=$4
    SNV_dir=$output
    INDEL_dir=$output
    CNV_dir=$output
    SV_dir=$output
    report_dir=$output/ANNOT/	
    mkdir -p $report_dir 
    cd $report_dir
    group=$sample

########################################################	
######		Reading run_info.txt and assigning to variables

    input=$( cat $run_info | grep -w '^INPUT_DIR' | cut -d '=' -f2)
    tool_info=$( cat $run_info | grep -w '^TOOL_INFO' | cut -d '=' -f2)
    
    type=$( cat $run_info | grep -w '^TYPE' | cut -d '=' -f2)
    script_path=$( cat $tool_info | grep -w '^WORKFLOW_PATH' | cut -d '=' -f2 )
    bedtools=$( cat $tool_info | grep -w '^BEDTOOLS' | cut -d '=' -f2 )
    master_gene_file=$( cat $tool_info | grep -w '^MASTER_GENE_FILE' | cut -d '=' -f2 )
    master_entrez_file=$( cat $tool_info | grep -w '^MASTER_ENTREZ_FILE' | cut -d '=' -f2 )
    email=$( cat $run_info | grep -w '^EMAIL' | cut -d '=' -f2)
    queue=$( cat $run_info | grep -w '^QUEUE' | cut -d '=' -f2)
    multi_sample=$( cat $run_info | grep -w '^MULTISAMPLE' | cut -d '=' -f2)
	Rsoft=$( cat $tool_info | grep -w '^R_SOFT' | cut -d '=' -f2 )
	somatic_calling=$( cat $tool_info | grep -w '^SOMATIC_CALLING' | cut -d '=' -f2 )
	export PATH=$Rsoft:$PATH
##############################################################		

		
    if [ $multi_sample != "YES" ]
    then
		echo "Single sample"
		cat $master_gene_file | cut -f4 > $report_dir/$sample.gene.temp
		if [ $type == exome ]
		then
			echo "Exome Analysis"
			### summarizing SNV files
			file=$SNV_dir/$sample.SNV.filtered.xls
			if [ ! -f $file ]
			then
				$script_path/errorlog.sh $file gene_summary.sh ERROR "not found"
				exit 1;
			fi	
			function=`awk -F '\t' '{ for(i=1;i<=NF;i++){ if ($i == "Effect") {print i} } }' $file`
			gene=`awk -F '\t' '{ for(i=1;i<=NF;i++){ if ($i == "geneList") {print i} } }' $file`
			
			cat $file | awk 'NR>2' | cut -f "$function","$gene" > $SNV_dir/$sample.SNV.tmp	
			for snv in SPLICE_SITE_ACCEPTOR SPLICE_SITE_DONOR START_LOST STOP_GAINED STOP_LOST RARE_AMINO_ACID NON_SYNONYMOUS_CODING SYNONYMOUS_START NON_SYNONYMOUS_START START_GAINED SYNONYMOUS_CODING SYNONYMOUS_STOP NON_SYNONYMOUS_STOP UTR_5_PRIME UTR_3_PRIME
			do
				cat $SNV_dir/$sample.SNV.tmp | grep "$snv" | cut -f2 | tr "," "\n" > $SNV_dir/in.$sample.$snv.tmp
				if [ ! -s $SNV_dir/in.$sample.$snv.tmp ]
				then
					echo "NOGENE"  >> $SNV_dir/in.$sample.$snv.tmp
				fi	
			done    
			Rscript $script_path/summary.SNV.r $report_dir/$sample.gene.temp $SNV_dir/in.$sample.SPLICE_SITE_ACCEPTOR.tmp $SNV_dir/in.$sample.SPLICE_SITE_DONOR.tmp $SNV_dir/in.$sample.START_LOST.tmp $SNV_dir/in.$sample.STOP_GAINED.tmp $SNV_dir/in.$sample.STOP_LOST.tmp $SNV_dir/in.$sample.RARE_AMINO_ACID.tmp $SNV_dir/in.$sample.NON_SYNONYMOUS_CODING.tmp $SNV_dir/in.$sample.SYNONYMOUS_START.tmp $SNV_dir/in.$sample.NON_SYNONYMOUS_START.tmp $SNV_dir/in.$sample.START_GAINED.tmp $SNV_dir/in.$sample.SYNONYMOUS_CODING.tmp $SNV_dir/in.$sample.SYNONYMOUS_STOP.tmp $SNV_dir/in.$sample.NON_SYNONYMOUS_STOP.tmp $SNV_dir/in.$sample.UTR_5_PRIME.tmp $SNV_dir/in.$sample.UTR_3_PRIME.tmp $SNV_dir/$sample.SPLICE_SITE_ACCEPTOR.tmp $SNV_dir/$sample.SPLICE_SITE_DONOR.tmp $SNV_dir/$sample.START_LOST.tmp $SNV_dir/$sample.STOP_GAINED.tmp $SNV_dir/$sample.STOP_LOST.tmp $SNV_dir/$sample.RARE_AMINO_ACID.tmp $SNV_dir/$sample.NON_SYNONYMOUS_CODING.tmp $SNV_dir/$sample.SYNONYMOUS_START.tmp $SNV_dir/$sample.NON_SYNONYMOUS_START.tmp $SNV_dir/$sample.START_GAINED.tmp $SNV_dir/$sample.SYNONYMOUS_CODING.tmp $SNV_dir/$sample.SYNONYMOUS_STOP.tmp $SNV_dir/$sample.NON_SYNONYMOUS_STOP.tmp $SNV_dir/$sample.UTR_5_PRIME.tmp $SNV_dir/$sample.UTR_3_PRIME.tmp
			
			join $SNV_dir/$sample.SPLICE_SITE_ACCEPTOR.tmp $SNV_dir/$sample.SPLICE_SITE_DONOR.tmp > $SNV_dir/$sample.join1.txt
			join $SNV_dir/$sample.join1.txt $SNV_dir/$sample.START_LOST.tmp > $SNV_dir/$sample.join2.txt
			join $SNV_dir/$sample.join2.txt $SNV_dir/$sample.STOP_GAINED.tmp > $SNV_dir/$sample.join3.txt
			join $SNV_dir/$sample.join3.txt $SNV_dir/$sample.STOP_LOST.tmp > $SNV_dir/$sample.join4.txt
			join $SNV_dir/$sample.join4.txt $SNV_dir/$sample.RARE_AMINO_ACID.tmp > $SNV_dir/$sample.join5.txt
			join $SNV_dir/$sample.join5.txt $SNV_dir/$sample.NON_SYNONYMOUS_CODING.tmp > $SNV_dir/$sample.join6.txt
			join $SNV_dir/$sample.join6.txt $SNV_dir/$sample.SYNONYMOUS_START.tmp > $SNV_dir/$sample.join7.txt
			join $SNV_dir/$sample.join7.txt $SNV_dir/$sample.NON_SYNONYMOUS_START.tmp > $SNV_dir/$sample.join8.txt
			join $SNV_dir/$sample.join8.txt $SNV_dir/$sample.START_GAINED.tmp > $SNV_dir/$sample.join9.txt
			join $SNV_dir/$sample.join9.txt $SNV_dir/$sample.SYNONYMOUS_CODING.tmp > $SNV_dir/$sample.join10.txt
			join $SNV_dir/$sample.join10.txt $SNV_dir/$sample.SYNONYMOUS_STOP.tmp > $SNV_dir/$sample.join11.txt
			join $SNV_dir/$sample.join11.txt $SNV_dir/$sample.NON_SYNONYMOUS_STOP.tmp > $SNV_dir/$sample.join12.txt
			join $SNV_dir/$sample.join12.txt $SNV_dir/$sample.UTR_5_PRIME.tmp > $SNV_dir/$sample.join13.txt
			join $SNV_dir/$sample.join13.txt $SNV_dir/$sample.UTR_3_PRIME.tmp > $SNV_dir/$sample.join14.txt
			cat $SNV_dir/$sample.join14.txt | tr " " "\t" > $SNV_dir/$sample.join15.txt
				
			touch $SNV_dir/$sample.SNV.summary
			echo -e "SPLICE_SITE_ACCEPTOR\tSPLICE_SITE_DONOR\tSTART_LOST\tSTOP_GAINED\tSTOP_LOST\tRARE_AMINO_ACID\tNON_SYNONYMOUS_CODING\tSYNONYMOUS_START\tNON_SYNONYMOUS_START\tSTART_GAINED\tSYNONYMOUS_CODING\tSYNONYMOUS_STOP\tNON_SYNONYMOUS_STOP\tUTR_5_PRIME\tUTR_3_PRIME" >> $SNV_dir/$sample.SNV.summary
			cat $SNV_dir/$sample.join15.txt >> $SNV_dir/$sample.SNV.summary
			Rscript $script_path/sum.cols.r $SNV_dir/$sample.SNV.summary $SNV_dir/$sample.SNV.sum	
			rm $SNV_dir/in.$sample.*.tmp $SNV_dir/$sample.*.tmp 
			rm $SNV_dir/$sample.join*.txt 
			#################################################################################################	
			### summarizing INDEL files
			file=$INDEL_dir/$sample.INDEL.filtered.xls
			if [ ! -f $file ]
			then
				$script_path/errorlog.sh $file gene_summary.sh ERROR "not found"
				exit 1;
			fi	
			function=`awk -F '\t' '{ for(i=1;i<=NF;i++){ if ($i == "Effect") {print i} } }' $file`
			gene=`awk -F '\t' '{ for(i=1;i<=NF;i++){ if ($i == "geneList") {print i} } }' $file`
			cat $file | awk 'NR>2' | cut -f "$function","$gene" > $INDEL_dir/$sample.INDEL.tmp

			for indel in EXON_DELETED FRAME_SHIFT CODON_CHANGE UTR_5_DELETED UTR_3_DELETED CODON_INSERTION CODON_CHANGE_PLUS_CODON_INSERTION CODON_DELETION CODON_CHANGE_PLUS_CODON_DELETION SPLICE_SITE_ACCEPTOR SPLICE_SITE_DONOR UTR_5_PRIME UTR_3_PRIME		
			do
				cat $INDEL_dir/$sample.INDEL.tmp | grep "$indel" | cut -f2 | tr "," "\n" > $INDEL_dir/in.$sample.$indel.tmp
				if [ ! -s $INDEL_dir/in.$sample.$indel.tmp ]
				then
					echo "NOGENE"  >> $INDEL_dir/in.$sample.$indel.tmp
				fi				
			done
			
			Rscript $script_path/summary.INDEL.r $report_dir/$sample.gene.temp $INDEL_dir/in.$sample.EXON_DELETED.tmp $INDEL_dir/in.$sample.FRAME_SHIFT.tmp $INDEL_dir/in.$sample.CODON_CHANGE.tmp $INDEL_dir/in.$sample.UTR_5_DELETED.tmp $INDEL_dir/in.$sample.UTR_3_DELETED.tmp $INDEL_dir/in.$sample.CODON_INSERTION.tmp $INDEL_dir/in.$sample.CODON_CHANGE_PLUS_CODON_INSERTION.tmp $INDEL_dir/in.$sample.CODON_DELETION.tmp $INDEL_dir/in.$sample.CODON_CHANGE_PLUS_CODON_DELETION.tmp $INDEL_dir/in.$sample.SPLICE_SITE_ACCEPTOR.tmp $INDEL_dir/in.$sample.SPLICE_SITE_DONOR.tmp $INDEL_dir/in.$sample.UTR_5_PRIME.tmp $INDEL_dir/in.$sample.UTR_3_PRIME.tmp $INDEL_dir/$sample.EXON_DELETED.tmp $INDEL_dir/$sample.FRAME_SHIFT.tmp $INDEL_dir/$sample.CODON_CHANGE.tmp $INDEL_dir/$sample.UTR_5_DELETED.tmp $INDEL_dir/$sample.UTR_3_DELETED.tmp $INDEL_dir/$sample.CODON_INSERTION.tmp $INDEL_dir/$sample.CODON_CHANGE_PLUS_CODON_INSERTION.tmp $INDEL_dir/$sample.CODON_DELETION.tmp $INDEL_dir/$sample.CODON_CHANGE_PLUS_CODON_DELETION.tmp $INDEL_dir/$sample.SPLICE_SITE_ACCEPTOR.tmp $INDEL_dir/$sample.SPLICE_SITE_DONOR.tmp $INDEL_dir/$sample.UTR_5_PRIME.tmp $INDEL_dir/$sample.UTR_3_PRIME.tmp
				
			join $INDEL_dir/$sample.EXON_DELETED.tmp $INDEL_dir/$sample.FRAME_SHIFT.tmp > $INDEL_dir/$sample.join1.txt
			join $INDEL_dir/$sample.join1.txt $INDEL_dir/$sample.CODON_CHANGE.tmp > $INDEL_dir/$sample.join2.txt
			join $INDEL_dir/$sample.join2.txt $INDEL_dir/$sample.UTR_5_DELETED.tmp > $INDEL_dir/$sample.join3.txt
			join $INDEL_dir/$sample.join3.txt $INDEL_dir/$sample.UTR_3_DELETED.tmp > $INDEL_dir/$sample.join4.txt
			join $INDEL_dir/$sample.join4.txt $INDEL_dir/$sample.CODON_INSERTION.tmp > $INDEL_dir/$sample.join5.txt
			join $INDEL_dir/$sample.join5.txt $INDEL_dir/$sample.CODON_CHANGE_PLUS_CODON_INSERTION.tmp > $INDEL_dir/$sample.join6.txt
			join $INDEL_dir/$sample.join6.txt $INDEL_dir/$sample.CODON_DELETION.tmp > $INDEL_dir/$sample.join7.txt
			join $INDEL_dir/$sample.join7.txt $INDEL_dir/$sample.CODON_CHANGE_PLUS_CODON_DELETION.tmp > $INDEL_dir/$sample.join8.txt
			join $INDEL_dir/$sample.join8.txt $INDEL_dir/$sample.SPLICE_SITE_ACCEPTOR.tmp > $INDEL_dir/$sample.join9.txt
			join $INDEL_dir/$sample.join9.txt $INDEL_dir/$sample.SPLICE_SITE_DONOR.tmp > $INDEL_dir/$sample.join10.txt
			join $INDEL_dir/$sample.join10.txt $INDEL_dir/$sample.UTR_5_PRIME.tmp > $INDEL_dir/$sample.join11.txt
			join $INDEL_dir/$sample.join11.txt $INDEL_dir/$sample.UTR_3_PRIME.tmp > $INDEL_dir/$sample.join12.txt
			cat $INDEL_dir/$sample.join12.txt | tr " " "\t" > $INDEL_dir/$sample.join13.txt
				
			touch $INDEL_dir/$sample.INDEL.summary
			echo -e "EXON_DELETED\tFRAME_SHIFT\tCODON_CHANGE\tUTR_5_DELETED\tUTR_3_DELETED\tCODON_INSERTION\tCODON_CHANGE_PLUS_CODON_INSERTION\tCODON_DELETION\tCODON_CHANGE_PLUS_CODON_DELETION\tSPLICE_SITE_ACCEPTOR\tSPLICE_SITE_DONOR\tUTR_5_PRIME\tUTR_3_PRIME" >> $INDEL_dir/$sample.INDEL.summary
			cat $INDEL_dir/$sample.join13.txt >> $INDEL_dir/$sample.INDEL.summary
			Rscript $script_path/sum.cols.r $INDEL_dir/$sample.INDEL.summary $INDEL_dir/$sample.INDEL.sum
			
			rm $INDEL_dir/in.$sample.*.tmp $INDEL_dir/$sample.*.tmp $INDEL_dir/$sample.join*.txt 

			#################################################################################################					### generating gene summary file
			cat $master_gene_file | awk '{print $4"\t"$1"\t"$2"\t"$3"\t"$5}' > $report_dir/$sample.GeneList.forsummary.txt
			cat $master_entrez_file | awk '{print $2}' > $report_dir/$sample.EntrezID.txt
			touch $report_dir/$sample.Gene.Summary.txt
			echo -e "\t\t\t\t\t\t\t\t\t\t$sample" >> $report_dir/$sample.Gene.Summary.txt
			echo -e "\t\t\t\t\t\t\t\tSNV_BREAKDOWN\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tINDEL_BREAKDOWN" >> $report_dir/$sample.Gene.Summary.txt
			echo -e "GENE\tCHROMOSOME\tSTART\tSTOP\tSTRAND\tENTREZ_GENE_ID\tTOTAL_SNVs\tTOTAL_INDELs\tSPLICE_SITE_ACCEPTOR\tSPLICE_SITE_DONOR\tSTART_LOST\tSTOP_GAINED\tSTOP_LOST\tRARE_AMINO_ACID\tNON_SYNONYMOUS_CODING\tSYNONYMOUS_START\tNON_SYNONYMOUS_START\tSTART_GAINED\tSYNONYMOUS_CODING\tSYNONYMOUS_STOP\tNON_SYNONYMOUS_STOP\tUTR_5_PRIME\tUTR_3_PRIME\tEXON_DELETED\tFRAME_SHIFT\tCODON_CHANGE\tUTR_5_DELETED\tUTR_3_DELETED\tCODON_INSERTION\tCODON_CHANGE_PLUS_CODON_INSERTION\tCODON_DELETION\tCODON_CHANGE_PLUS_CODON_DELETION\tSPLICE_SITE_ACCEPTOR\tSPLICE_SITE_DONOR\tUTR_5_PRIME\tUTR_3_PRIME" >> $report_dir/$sample.Gene.Summary.txt
			sed -i '1d' $SNV_dir/$sample.SNV.summary
			sed -i '1d' $INDEL_dir/$sample.INDEL.summary
			cat $SNV_dir/$sample.SNV.summary | cut -f2-16 > $SNV_dir/$sample.SNV.tmp
			cat $INDEL_dir/$sample.INDEL.summary | cut -f2-14 > $INDEL_dir/$sample.INDEL.tmp
			
			paste $report_dir/$sample.GeneList.forsummary.txt $report_dir/$sample.EntrezID.txt $SNV_dir/$sample.SNV.sum $INDEL_dir/$sample.INDEL.sum $SNV_dir/$sample.SNV.tmp $INDEL_dir/$sample.INDEL.tmp >> $report_dir/$sample.Gene.Summary.txt
			
			rm $SNV_dir/$sample.SNV.tmp $INDEL_dir/$sample.INDEL.tmp $report_dir/$sample.GeneList.forsummary.txt $INDEL_dir/$sample.INDEL.summary $INDEL_dir/$sample.INDEL.sum $SNV_dir/$sample.SNV.summary $SNV_dir/$sample.SNV.sum $report_dir/$sample.EntrezID.txt
			rm $report_dir/$sample.gene.temp
				
			#################################################################################################	
		else
			echo "Whole Genome Anlaysis"
			### summarizing SNV files
			file=$SNV_dir/$sample.SNV.filtered.xls
			if [ ! -f $file ]
			then
				$script_path/errorlog.sh $file gene_summary.sh ERROR "not found"
				exit 1;
			fi	
			function=`awk -F '\t' '{ for(i=1;i<=NF;i++){ if ($i == "Effect") {print i} } }' $file`
			gene=`awk -F '\t' '{ for(i=1;i<=NF;i++){ if ($i == "geneList") {print i} } }' $file`
			cat $file | awk 'NR>2' | cut -f "$function","$gene" > $SNV_dir/$sample.SNV.tmp	
			for snv in SPLICE_SITE_ACCEPTOR SPLICE_SITE_DONOR START_LOST STOP_GAINED STOP_LOST RARE_AMINO_ACID NON_SYNONYMOUS_CODING SYNONYMOUS_START NON_SYNONYMOUS_START START_GAINED SYNONYMOUS_CODING SYNONYMOUS_STOP NON_SYNONYMOUS_STOP UTR_5_PRIME UTR_3_PRIME
			do
				cat $SNV_dir/$sample.SNV.tmp | grep "$snv" | cut -f2 | tr "," "\n" > $SNV_dir/in.$sample.$snv.tmp
				if [ ! -s $SNV_dir/in.$sample.$snv.tmp ]
				then
					echo "NOGENE"  >> $SNV_dir/in.$sample.$snv.tmp
				fi	
			done    
			
			Rscript $script_path/summary.SNV.r $report_dir/$sample.gene.temp $SNV_dir/in.$sample.SPLICE_SITE_ACCEPTOR.tmp $SNV_dir/in.$sample.SPLICE_SITE_DONOR.tmp $SNV_dir/in.$sample.START_LOST.tmp $SNV_dir/in.$sample.STOP_GAINED.tmp $SNV_dir/in.$sample.STOP_LOST.tmp $SNV_dir/in.$sample.RARE_AMINO_ACID.tmp $SNV_dir/in.$sample.NON_SYNONYMOUS_CODING.tmp $SNV_dir/in.$sample.SYNONYMOUS_START.tmp $SNV_dir/in.$sample.NON_SYNONYMOUS_START.tmp $SNV_dir/in.$sample.START_GAINED.tmp $SNV_dir/in.$sample.SYNONYMOUS_CODING.tmp $SNV_dir/in.$sample.SYNONYMOUS_STOP.tmp $SNV_dir/in.$sample.NON_SYNONYMOUS_STOP.tmp $SNV_dir/in.$sample.UTR_5_PRIME.tmp $SNV_dir/in.$sample.UTR_3_PRIME.tmp $SNV_dir/$sample.SPLICE_SITE_ACCEPTOR.tmp $SNV_dir/$sample.SPLICE_SITE_DONOR.tmp $SNV_dir/$sample.START_LOST.tmp $SNV_dir/$sample.STOP_GAINED.tmp $SNV_dir/$sample.STOP_LOST.tmp $SNV_dir/$sample.RARE_AMINO_ACID.tmp $SNV_dir/$sample.NON_SYNONYMOUS_CODING.tmp $SNV_dir/$sample.SYNONYMOUS_START.tmp $SNV_dir/$sample.NON_SYNONYMOUS_START.tmp $SNV_dir/$sample.START_GAINED.tmp $SNV_dir/$sample.SYNONYMOUS_CODING.tmp $SNV_dir/$sample.SYNONYMOUS_STOP.tmp $SNV_dir/$sample.NON_SYNONYMOUS_STOP.tmp $SNV_dir/$sample.UTR_5_PRIME.tmp $SNV_dir/$sample.UTR_3_PRIME.tmp
			
			join $SNV_dir/$sample.SPLICE_SITE_ACCEPTOR.tmp $SNV_dir/$sample.SPLICE_SITE_DONOR.tmp > $SNV_dir/$sample.join1.txt
			join $SNV_dir/$sample.join1.txt $SNV_dir/$sample.START_LOST.tmp > $SNV_dir/$sample.join2.txt
			join $SNV_dir/$sample.join2.txt $SNV_dir/$sample.STOP_GAINED.tmp > $SNV_dir/$sample.join3.txt
			join $SNV_dir/$sample.join3.txt $SNV_dir/$sample.STOP_LOST.tmp > $SNV_dir/$sample.join4.txt
			join $SNV_dir/$sample.join4.txt $SNV_dir/$sample.RARE_AMINO_ACID.tmp > $SNV_dir/$sample.join5.txt
			join $SNV_dir/$sample.join5.txt $SNV_dir/$sample.NON_SYNONYMOUS_CODING.tmp > $SNV_dir/$sample.join6.txt
			join $SNV_dir/$sample.join6.txt $SNV_dir/$sample.SYNONYMOUS_START.tmp > $SNV_dir/$sample.join7.txt
			join $SNV_dir/$sample.join7.txt $SNV_dir/$sample.NON_SYNONYMOUS_START.tmp > $SNV_dir/$sample.join8.txt
			join $SNV_dir/$sample.join8.txt $SNV_dir/$sample.START_GAINED.tmp > $SNV_dir/$sample.join9.txt
			join $SNV_dir/$sample.join9.txt $SNV_dir/$sample.SYNONYMOUS_CODING.tmp > $SNV_dir/$sample.join10.txt
			join $SNV_dir/$sample.join10.txt $SNV_dir/$sample.SYNONYMOUS_STOP.tmp > $SNV_dir/$sample.join11.txt
			join $SNV_dir/$sample.join11.txt $SNV_dir/$sample.NON_SYNONYMOUS_STOP.tmp > $SNV_dir/$sample.join12.txt
			join $SNV_dir/$sample.join12.txt $SNV_dir/$sample.UTR_5_PRIME.tmp > $SNV_dir/$sample.join13.txt
			join $SNV_dir/$sample.join13.txt $SNV_dir/$sample.UTR_3_PRIME.tmp > $SNV_dir/$sample.join14.txt
			cat $SNV_dir/$sample.join14.txt | tr " " "\t" > $SNV_dir/$sample.join15.txt
				
			touch $SNV_dir/$sample.SNV.summary
			echo -e "SPLICE_SITE_ACCEPTOR\tSPLICE_SITE_DONOR\tSTART_LOST\tSTOP_GAINED\tSTOP_LOST\tRARE_AMINO_ACID\tNON_SYNONYMOUS_CODING\tSYNONYMOUS_START\tNON_SYNONYMOUS_START\tSTART_GAINED\tSYNONYMOUS_CODING\tSYNONYMOUS_STOP\tNON_SYNONYMOUS_STOP\tUTR_5_PRIME\tUTR_3_PRIME" >> $SNV_dir/$sample.SNV.summary
			cat $SNV_dir/$sample.join15.txt >> $SNV_dir/$sample.SNV.summary
			Rscript $script_path/sum.cols.r $SNV_dir/$sample.SNV.summary $SNV_dir/$sample.SNV.sum
				
			rm $SNV_dir/$sample.*.tmp $SNV_dir/in.$sample.*.tmp $SNV_dir/$sample.join*.txt 

			#################################################################################################	
				### summarizing INDEL files
			file=$INDEL_dir/$sample.INDEL.filtered.xls
			if [ ! -f $file ]
			then
				$script_path/errorlog.sh $file gene_summary.sh ERROR "not found"
				exit 1;
			fi	
			function=`awk -F '\t' '{ for(i=1;i<=NF;i++){ if ($i == "Effect") {print i} } }' $file`
			gene=`awk -F '\t' '{ for(i=1;i<=NF;i++){ if ($i == "geneList") {print i} } }' $file`
			cat $file | awk 'NR>2' | cut -f "$function","$gene" > $INDEL_dir/$sample.INDEL.tmp

			for indel in EXON_DELETED FRAME_SHIFT CODON_CHANGE UTR_5_DELETED UTR_3_DELETED CODON_INSERTION CODON_CHANGE_PLUS_CODON_INSERTION CODON_DELETION CODON_CHANGE_PLUS_CODON_DELETION SPLICE_SITE_ACCEPTOR SPLICE_SITE_DONOR UTR_5_PRIME UTR_3_PRIME		
			do
				cat $INDEL_dir/$sample.INDEL.tmp | grep "$indel" | cut -f2 | tr "," "\n" > $INDEL_dir/in.$sample.$indel.tmp
				if [ ! -s $INDEL_dir/$sample.$indel.tmp ]
				then
					echo "NOGENE"  >> $INDEL_dir/in.$sample.$indel.tmp
				fi				
			done
				
			Rscript $script_path/summary.INDEL.r $report_dir/$sample.gene.temp $INDEL_dir/in.$sample.EXON_DELETED.tmp $INDEL_dir/in.$sample.FRAME_SHIFT.tmp $INDEL_dir/in.$sample.CODON_CHANGE.tmp $INDEL_dir/in.$sample.UTR_5_DELETED.tmp $INDEL_dir/in.$sample.UTR_3_DELETED.tmp $INDEL_dir/in.$sample.CODON_INSERTION.tmp $INDEL_dir/in.$sample.CODON_CHANGE_PLUS_CODON_INSERTION.tmp $INDEL_dir/in.$sample.CODON_DELETION.tmp $INDEL_dir/in.$sample.CODON_CHANGE_PLUS_CODON_DELETION.tmp $INDEL_dir/in.$sample.SPLICE_SITE_ACCEPTOR.tmp $INDEL_dir/in.$sample.SPLICE_SITE_DONOR.tmp $INDEL_dir/in.$sample.UTR_5_PRIME.tmp $INDEL_dir/in.$sample.UTR_3_PRIME.tmp $INDEL_dir/$sample.EXON_DELETED.tmp $INDEL_dir/$sample.FRAME_SHIFT.tmp $INDEL_dir/$sample.CODON_CHANGE.tmp $INDEL_dir/$sample.UTR_5_DELETED.tmp $INDEL_dir/$sample.UTR_3_DELETED.tmp $INDEL_dir/$sample.CODON_INSERTION.tmp $INDEL_dir/$sample.CODON_CHANGE_PLUS_CODON_INSERTION.tmp $INDEL_dir/$sample.CODON_DELETION.tmp $INDEL_dir/$sample.CODON_CHANGE_PLUS_CODON_DELETION.tmp $INDEL_dir/$sample.SPLICE_SITE_ACCEPTOR.tmp $INDEL_dir/$sample.SPLICE_SITE_DONOR.tmp $INDEL_dir/$sample.UTR_5_PRIME.tmp $INDEL_dir/$sample.UTR_3_PRIME.tmp
			
			join $INDEL_dir/$sample.EXON_DELETED.tmp $INDEL_dir/$sample.FRAME_SHIFT.tmp > $INDEL_dir/$sample.join1.txt
			join $INDEL_dir/$sample.join1.txt $INDEL_dir/$sample.CODON_CHANGE.tmp > $INDEL_dir/$sample.join2.txt
			join $INDEL_dir/$sample.join2.txt $INDEL_dir/$sample.UTR_5_DELETED.tmp > $INDEL_dir/$sample.join3.txt
			join $INDEL_dir/$sample.join3.txt $INDEL_dir/$sample.UTR_3_DELETED.tmp > $INDEL_dir/$sample.join4.txt
			join $INDEL_dir/$sample.join4.txt $INDEL_dir/$sample.CODON_INSERTION.tmp > $INDEL_dir/$sample.join5.txt
			join $INDEL_dir/$sample.join5.txt $INDEL_dir/$sample.CODON_CHANGE_PLUS_CODON_INSERTION.tmp > $INDEL_dir/$sample.join6.txt
			join $INDEL_dir/$sample.join6.txt $INDEL_dir/$sample.CODON_DELETION.tmp > $INDEL_dir/$sample.join7.txt
			join $INDEL_dir/$sample.join7.txt $INDEL_dir/$sample.CODON_CHANGE_PLUS_CODON_DELETION.tmp > $INDEL_dir/$sample.join8.txt
			join $INDEL_dir/$sample.join8.txt $INDEL_dir/$sample.SPLICE_SITE_ACCEPTOR.tmp > $INDEL_dir/$sample.join9.txt
			join $INDEL_dir/$sample.join9.txt $INDEL_dir/$sample.SPLICE_SITE_DONOR.tmp > $INDEL_dir/$sample.join10.txt
			join $INDEL_dir/$sample.join10.txt $INDEL_dir/$sample.UTR_5_PRIME.tmp > $INDEL_dir/$sample.join11.txt
			join $INDEL_dir/$sample.join11.txt $INDEL_dir/$sample.UTR_3_PRIME.tmp > $INDEL_dir/$sample.join12.txt

			cat $INDEL_dir/$sample.join12.txt | tr " " "\t" > $INDEL_dir/$sample.join13.txt
			
			touch $INDEL_dir/$sample.INDEL.summary
			echo -e "EXON_DELETED\tFRAME_SHIFT\tCODON_CHANGE\tUTR_5_DELETED\tUTR_3_DELETED\tCODON_INSERTION\tCODON_CHANGE_PLUS_CODON_INSERTION\tCODON_DELETION\tCODON_CHANGE_PLUS_CODON_DELETION\tSPLICE_SITE_ACCEPTOR\tSPLICE_SITE_DONOR\tUTR_5_PRIME\tUTR_3_PRIME" >> $INDEL_dir/$sample.INDEL.summary
			cat $INDEL_dir/$sample.join13.txt >> $INDEL_dir/$sample.INDEL.summary
			Rscript $script_path/sum.cols.r $INDEL_dir/$sample.INDEL.summary $INDEL_dir/$sample.INDEL.sum
			
			rm $INDEL_dir/$sample.*.tmp $INDEL_dir/in.$sample.*.tmp $INDEL_dir/$sample.join*.txt 
				#################################################################################################	
				### summarizing CNV files
			file=$INDEL_dir/$sample.INDEL.filtered.xls
			if [ ! -f $file ]
			then
				$script_path/errorlog.sh $file gene_summary.sh ERROR "not found"
				exit 1;
			fi	
			function=`awk -F '\t' '{ for(i=1;i<=NF;i++){ if ($i == "Effect") {print i} } }' $file`
			gene=`awk -F '\t' '{ for(i=1;i<=NF;i++){ if ($i == "geneList") {print i} } }' $file`
			cat $file | awk 'NR>2' | cut -f "$function","$gene" > $SNV_dir/$sample.SNV.tmp	

			cat $CNV_dir/ANNOT/$sample.CNV.annotated.txt | tr "_" "\t" | cut -f4,10 | grep DUP > $CNV_dir/$sample.DUP.tmp
			cat $CNV_dir/ANNOT/$sample.CNV.annotated.txt | tr "_" "\t" | cut -f4,10 | grep DEL > $CNV_dir/$sample.DEL.tmp
			
			if [ ! -s $CNV_dir/$sample.DUP.tmp ]
			then
				echo "NOGENE" >> $CNV_dir/$sample.DUP.tmp
			fi
			if [ ! -s $CNV_dir/$sample.DEL.tmp ]
			then
				echo "NOGENE" >> $CNV_dir/$sample.DEL.tmp
			fi

			Rscript $script_path/summary.CNV.r $report_dir/$sample.gene.temp $CNV_dir/$sample.DEL.tmp $CNV_dir/$sample.DUP.tmp $CNV_dir/$sample.DEL.txt $CNV_dir/$sample.DUP.txt
			
			join $CNV_dir/$sample.DEL.txt $CNV_dir/$sample.DUP.txt > $CNV_dir/$sample.join.txt
			cat $CNV_dir/$sample.join.txt | tr " " "\t" > $CNV_dir/$sample.join1.txt
			
			touch $CNV_dir/$sample.CNV.summary
			echo -e "DEL\tDUP" >> $CNV_dir/$sample.CNV.summary
			cat $CNV_dir/$sample.join1.txt >> $CNV_dir/$sample.CNV.summary
			Rscript $script_path/sum.cols.r $CNV_dir/$sample.CNV.summary $CNV_dir/$sample.CNV.sum
			
			rm $CNV_dir/$sample.*.tmp $CNV_dir/$sample.join*.txt $CNV_dir/$sample.DEL.txt $CNV_dir/$sample.DUP.txt
	#################################################################################################					### summarizing SV files
			cat $SV_dir/ANNOT/$sample.SV.annotated.txt | sed -e '/NA_/s//NOGENE_/g' -e '/_NA/s//_NOGENE/g' | grep INV | tr "_" "\t" | cut -f9,10,11 | tr " " "\t" > $SV_dir/$sample.INV.tmp
			cat $SV_dir/ANNOT/$sample.SV.annotated.txt | sed -e '/NA_/s//NOGENE_/g' -e '/_NA/s//_NOGENE/g' | grep INS | tr "_" "\t" | cut -f9,10,11 | tr " " "\t" > $SV_dir/$sample.INS.tmp
			cat $SV_dir/ANNOT/$sample.SV.annotated.txt | sed -e '/NA_/s//NOGENE_/g' -e '/_NA/s//_NOGENE/g' | grep DEL | tr "_" "\t" | cut -f9,10,11 | tr " " "\t" > $SV_dir/$sample.DEL.tmp
			cat $SV_dir/ANNOT/$sample.SV.annotated.txt | sed -e '/NA_/s//NOGENE_/g' -e '/_NA/s//_NOGENE/g' | grep ITX | tr "_" "\t" | cut -f9,10,11 | tr " " "\t" > $SV_dir/$sample.ITX.tmp
			cat $SV_dir/ANNOT/$sample.SV.annotated.txt | sed -e '/NA_/s//NOGENE_/g' -e '/_NA/s//_NOGENE/g' | grep CTX | tr "_" "\t" | cut -f9,10,11 | tr " " "\t" > $SV_dir/$sample.CTX.tmp
			
			if [ ! -s $SV_dir/$sample.INV.tmp ]
			then
				echo "NOGENE" >> $SV_dir/$sample.INV.tmp
			fi
			if [ ! -s $SV_dir/$sample.INS.tmp ]
			then
				echo "NOGENE" >> $SV_dir/$sample.INS.tmp
			fi
			if [ ! -s $SV_dir/$sample.DEL.tmp ]
			then
				echo "NOGENE" >> $SV_dir/$sample.DEL.tmp
			fi
			if [ ! -s $SV_dir/$sample.ITX.tmp ]
			then
				echo "NOGENE" >> $SV_dir/$sample.ITX.tmp
			fi
			if [ ! -s $SV_dir/$sample.CTX.tmp ]
			then
				echo "NOGENE" >> $SV_dir/$sample.CTX.tmp
			fi

			Rscript $script_path/summary.SV.r $report_dir/$sample.gene.temp $SV_dir/$sample.ITX.tmp $SV_dir/$sample.INV.tmp $SV_dir/$sample.DEL.tmp $SV_dir/$sample.INS.tmp $SV_dir/$sample.CTX.tmp $SV_dir/$sample.ITX.txt $SV_dir/$sample.INV.txt $SV_dir/$sample.DEL.txt $SV_dir/$sample.INS.txt $SV_dir/$sample.CTX.txt
			
			join $SV_dir/$sample.ITX.txt $SV_dir/$sample.INV.txt > $SV_dir/$sample.join1.txt
			join $SV_dir/$sample.join1.txt $SV_dir/$sample.DEL.txt > $SV_dir/$sample.join2.txt
			join $SV_dir/$sample.join2.txt $SV_dir/$sample.INS.txt > $SV_dir/$sample.join3.txt
			join $SV_dir/$sample.join3.txt $SV_dir/$sample.CTX.txt > $SV_dir/$sample.join4.txt
			cat $SV_dir/$sample.join4.txt | tr " " "\t" > $SV_dir/$sample.join5.txt
			
			touch $SV_dir/$sample.SV.summary
			echo -e "ITX\tINV\tDEL\tINS\tCTX" >> $SV_dir/$sample.SV.summary
			cat $SV_dir/$sample.join5.txt >> $SV_dir/$sample.SV.summary
			Rscript $script_path/sum.cols.r $SV_dir/$sample.SV.summary $SV_dir/$sample.SV.sum
			
			rm $SV_dir/$sample.*.tmp $SV_dir/$sample.join*.txt $SV_dir/$sample.ITX.txt $SV_dir/$sample.INV.txt $SV_dir/$sample.DEL.txt $SV_dir/$sample.INS.txt $SV_dir/$sample.CTX.txt 
	#################################################################################################					### generating gene summary file
			cat $master_gene_file | awk '{print $4"\t"$1"\t"$2"\t"$3"\t"$5}' > $report_dir/$sample.GeneList.forsummary.txt
			cat $master_entrez_file | awk '{print $2}' > $report_dir/$sample.EntrezID.txt
			touch $report_dir/$sample.Gene.Summary.txt
			echo -e "\t\t\t\t\t\t\t\t\t\t$sample" >> $report_dir/$sample.Gene.Summary.txt
			echo -e "\t\t\t\t\t\t\t\t\t\tSNV_BREAKDOWN\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tINDEL_BREAKDOWN\t\t\t\t\t\t\t\t\t\t\t\t\tCNV_BREAKDOWN\t\tSV_BREAKDOWN\t\t\t\t" >> $report_dir/$sample.Gene.Summary.txt
			echo -e "GENE\tCHROMOSOME\tSTART\tSTOP\tSTRAND\tENTREZ_GENE_ID\tTOTAL_SNVs\tTOTAL_INDELs\tTOTAL_CNVs\tTOTAL_SVs\tSPLICE_SITE_ACCEPTOR\tSPLICE_SITE_DONOR\tSTART_LOST\tSTOP_GAINED\tSTOP_LOST\tRARE_AMINO_ACID\tNON_SYNONYMOUS_CODING\tSYNONYMOUS_START\tNON_SYNONYMOUS_START\tSTART_GAINED\tSYNONYMOUS_CODING\tSYNONYMOUS_STOP\tNON_SYNONYMOUS_STOP\tUTR_5_PRIME\tUTR_3_PRIME\tEXON_DELETED\tFRAME_SHIFT\tCODON_CHANGE\tUTR_5_DELETED\tUTR_3_DELETED\tCODON_INSERTION\tCODON_CHANGE_PLUS_CODON_INSERTION\tCODON_DELETION\tCODON_CHANGE_PLUS_CODON_DELETION\tSPLICE_SITE_ACCEPTOR\tSPLICE_SITE_DONOR\tUTR_5_PRIME\tUTR_3_PRIME\tDELETION\tDUPLICATION\tITX\tINV\tDEL\tINS\tCTX" >> $report_dir/$sample.Gene.Summary.txt

			sed -i '1d' $SNV_dir/$sample.SNV.summary
			sed -i '1d' $INDEL_dir/$sample.INDEL.summary
			sed -i '1d' $CNV_dir/$sample.CNV.summary
			sed -i '1d' $SV_dir/$sample.SV.summary
			cat $SNV_dir/$sample.SNV.summary | cut -f2-16 > $SNV_dir/$sample.SNV.tmp
			cat $INDEL_dir/$sample.INDEL.summary | cut -f2-14 > $INDEL_dir/$sample.INDEL.tmp
			cat $CNV_dir/$sample.CNV.summary | cut -f2,3 > $CNV_dir/$sample.CNV.tmp
			cat $SV_dir/$sample.SV.summary | cut -f2-6 > $SV_dir/$sample.SV.tmp
			
			paste $report_dir/$sample.GeneList.forsummary.txt $report_dir/$sample.EntrezID.txt $SNV_dir/$sample.SNV.sum $INDEL_dir/$sample.INDEL.sum $CNV_dir/$sample.CNV.sum $SV_dir/$sample.SV.sum $SNV_dir/$sample.SNV.tmp $INDEL_dir/$sample.INDEL.tmp $CNV_dir/$sample.CNV.tmp $SV_dir/$sample.SV.tmp >> $report_dir/$sample.Gene.Summary.txt
			
			rm $SNV_dir/$sample.SNV.tmp $INDEL_dir/$sample.INDEL.tmp $CNV_dir/$sample.CNV.tmp $SV_dir/$sample.SV.tmp $report_dir/$sample.GeneList.forsummary.txt $SV_dir/$sample.SV.summary 
			rm $SV_dir/$sample.SV.sum $CNV_dir/$sample.CNV.summary $CNV_dir/$sample.CNV.sum $INDEL_dir/$sample.INDEL.summary $INDEL_dir/$sample.INDEL.sum $SNV_dir/$sample.SNV.summary $SNV_dir/$sample.SNV.sum $report_dir/$sample.EntrezID.txt
			rm $report_dir/$sample.gene.temp
		fi	
	else
		echo "Multi sample"
		if [ $type == "exome" ]
		then
			sample_info=$( cat $run_info | grep -w '^SAMPLE_INFO' | cut -d '=' -f2 )
			samples=$( cat $sample_info | grep -w "^$group" | cut -d '=' -f2 )
			let num_tumor=`echo $samples|tr " " "\n"|wc -l`-1
			tumor_list=`echo $samples | tr " " "\n" | tail -$num_tumor`
			if [ $somatic_calling == "NO" ]
			then
				tumor_list=$samples
			fi	
			for sample in $tumor_list    
			do
				if [ $somatic_calling == "YES" ]
				then
					sampe=$group.$sample
				else
					sampe=$sample
				fi	
				cat $master_gene_file | cut -f4 > $report_dir/$sampe.gene.temp
				if [ $somatic_calling == "YES" ]
				then
					file=$SNV_dir/TUMOR.$group.SNV.filtered.xls
				else
					file=$SNV_dir/$group.SNV.filtered.xls
				fi
				
				if [ ! -f $file ]
				then
					$script_path/errorlog.sh $file gene_summary.sh ERROR "not found"
					exit 1;
				fi	
				function=`awk -F '\t' '{ for(i=1;i<=NF;i++){ if ($i == "Effect") {print i} } }' $file`
				gene=`awk -F '\t' '{ for(i=1;i<=NF;i++){ if ($i == "geneList") {print i} } }' $file`
                sam=`awk -v sam=$sample -F '\t' '{ for(i=1;i<=NF;i++){ if ($i == sam) {print i} } }' $file`
				
				cat $file | awk 'NR>2' | awk -v sam=$sam '$sam !~ /n\/a/'| cut -f "$function","$gene" > $SNV_dir/$sampe.SNV.tmp
				for snv in SPLICE_SITE_ACCEPTOR SPLICE_SITE_DONOR START_LOST STOP_GAINED STOP_LOST RARE_AMINO_ACID NON_SYNONYMOUS_CODING SYNONYMOUS_START NON_SYNONYMOUS_START START_GAINED SYNONYMOUS_CODING SYNONYMOUS_STOP NON_SYNONYMOUS_STOP UTR_5_PRIME UTR_3_PRIME
				do
					cat $SNV_dir/$sampe.SNV.tmp | grep "$snv" | cut -f2 | tr "," "\n" > $SNV_dir/in.$sampe.$snv.tmp
					if [ ! -s $SNV_dir/$sampe.$snv.tmp ]
					then
						echo "NOGENE"  >> $SNV_dir/in.$sampe.$snv.tmp
					fi	
				done    
			
				Rscript $script_path/summary.SNV.r $report_dir/$sampe.gene.temp $SNV_dir/in.$sampe.SPLICE_SITE_ACCEPTOR.tmp $SNV_dir/in.$sampe.SPLICE_SITE_DONOR.tmp $SNV_dir/in.$sampe.START_LOST.tmp $SNV_dir/in.$sampe.STOP_GAINED.tmp $SNV_dir/in.$sampe.STOP_LOST.tmp $SNV_dir/in.$sampe.RARE_AMINO_ACID.tmp $SNV_dir/in.$sampe.NON_SYNONYMOUS_CODING.tmp $SNV_dir/in.$sampe.SYNONYMOUS_START.tmp $SNV_dir/in.$sampe.NON_SYNONYMOUS_START.tmp $SNV_dir/in.$sampe.START_GAINED.tmp $SNV_dir/in.$sampe.SYNONYMOUS_CODING.tmp $SNV_dir/in.$sampe.SYNONYMOUS_STOP.tmp $SNV_dir/in.$sampe.NON_SYNONYMOUS_STOP.tmp $SNV_dir/in.$sampe.UTR_5_PRIME.tmp $SNV_dir/in.$sampe.UTR_3_PRIME.tmp $SNV_dir/$sampe.SPLICE_SITE_ACCEPTOR.tmp $SNV_dir/$sampe.SPLICE_SITE_DONOR.tmp $SNV_dir/$sampe.START_LOST.tmp $SNV_dir/$sampe.STOP_GAINED.tmp $SNV_dir/$sampe.STOP_LOST.tmp $SNV_dir/$sampe.RARE_AMINO_ACID.tmp $SNV_dir/$sampe.NON_SYNONYMOUS_CODING.tmp $SNV_dir/$sampe.SYNONYMOUS_START.tmp $SNV_dir/$sampe.NON_SYNONYMOUS_START.tmp $SNV_dir/$sampe.START_GAINED.tmp $SNV_dir/$sampe.SYNONYMOUS_CODING.tmp $SNV_dir/$sampe.SYNONYMOUS_STOP.tmp $SNV_dir/$sampe.NON_SYNONYMOUS_STOP.tmp $SNV_dir/$sampe.UTR_5_PRIME.tmp $SNV_dir/$sampe.UTR_3_PRIME.tmp
			
				join $SNV_dir/$sampe.SPLICE_SITE_ACCEPTOR.tmp $SNV_dir/$sampe.SPLICE_SITE_DONOR.tmp > $SNV_dir/$sampe.join1.txt
				join $SNV_dir/$sampe.join1.txt $SNV_dir/$sampe.START_LOST.tmp > $SNV_dir/$sampe.join2.txt
				join $SNV_dir/$sampe.join2.txt $SNV_dir/$sampe.STOP_GAINED.tmp > $SNV_dir/$sampe.join3.txt
				join $SNV_dir/$sampe.join3.txt $SNV_dir/$sampe.STOP_LOST.tmp > $SNV_dir/$sampe.join4.txt
				join $SNV_dir/$sampe.join4.txt $SNV_dir/$sampe.RARE_AMINO_ACID.tmp > $SNV_dir/$sampe.join5.txt
				join $SNV_dir/$sampe.join5.txt $SNV_dir/$sampe.NON_SYNONYMOUS_CODING.tmp > $SNV_dir/$sampe.join6.txt
				join $SNV_dir/$sampe.join6.txt $SNV_dir/$sampe.SYNONYMOUS_START.tmp > $SNV_dir/$sampe.join7.txt
				join $SNV_dir/$sampe.join7.txt $SNV_dir/$sampe.NON_SYNONYMOUS_START.tmp > $SNV_dir/$sampe.join8.txt
				join $SNV_dir/$sampe.join8.txt $SNV_dir/$sampe.START_GAINED.tmp > $SNV_dir/$sampe.join9.txt
				join $SNV_dir/$sampe.join9.txt $SNV_dir/$sampe.SYNONYMOUS_CODING.tmp > $SNV_dir/$sampe.join10.txt
				join $SNV_dir/$sampe.join10.txt $SNV_dir/$sampe.SYNONYMOUS_STOP.tmp > $SNV_dir/$sampe.join11.txt
				join $SNV_dir/$sampe.join11.txt $SNV_dir/$sampe.NON_SYNONYMOUS_STOP.tmp > $SNV_dir/$sampe.join12.txt
				join $SNV_dir/$sampe.join12.txt $SNV_dir/$sampe.UTR_5_PRIME.tmp > $SNV_dir/$sampe.join13.txt
				join $SNV_dir/$sampe.join13.txt $SNV_dir/$sampe.UTR_3_PRIME.tmp > $SNV_dir/$sampe.join14.txt
				cat $SNV_dir/$sampe.join14.txt | tr " " "\t" > $SNV_dir/$sampe.join15.txt
					
				touch $SNV_dir/$sampe.SNV.summary
				echo -e "\tSPLICE_SITE_ACCEPTOR\tSPLICE_SITE_DONOR\tSTART_LOST\tSTOP_GAINED\tSTOP_LOST\tRARE_AMINO_ACID\tNON_SYNONYMOUS_CODING\tSYNONYMOUS_START\tNON_SYNONYMOUS_START\tSTART_GAINED\tSYNONYMOUS_CODING\tSYNONYMOUS_STOP\tNON_SYNONYMOUS_STOP\tUTR_5_PRIME\tUTR_3_PRIME" >> $SNV_dir/$sampe.SNV.summary
				cat $SNV_dir/$sampe.join15.txt >> $SNV_dir/$sampe.SNV.summary
				Rscript $script_path/sum.cols.r $SNV_dir/$sampe.SNV.summary $SNV_dir/$sampe.SNV.sum
					
				rm $SNV_dir/$sampe.*.tmp $SNV_dir/in.$sampe.*.tmp $SNV_dir/$sampe.join*.txt

				#################################################################################################	
				### summarizing INDEL files           
				if [ $somatic_calling == "YES" ]
				then
					file=$INDEL_dir/TUMOR.$group.INDEL.filtered.xls
				else
					file=$INDEL_dir/$group.INDEL.filtered.xls
				fi
				if [ ! -f $file ]
				then
					$script_path/errorlog.sh $file gene_summary.sh ERROR "not found"
					exit 1;
				fi	
				function=`awk -F '\t' '{ for(i=1;i<=NF;i++){ if ($i == "Effect") {print i} } }' $file`
				gene=`awk -F '\t' '{ for(i=1;i<=NF;i++){ if ($i == "geneList") {print i} } }' $file`
                                sam=`awk -v sam=$sample -F '\t' '{ for(i=1;i<=NF;i++){ if ($i == sam) {print i} } }' $file`
				cat $file | awk 'NR>2' | awk -v sam=$sam '$sam !~ /n\/a/'|  cut -f "$function","$gene" > $INDEL_dir/$sampe.INDEL.tmp

				for indel in EXON_DELETED FRAME_SHIFT CODON_CHANGE UTR_5_DELETED UTR_3_DELETED CODON_INSERTION CODON_CHANGE_PLUS_CODON_INSERTION CODON_DELETION CODON_CHANGE_PLUS_CODON_DELETION SPLICE_SITE_ACCEPTOR SPLICE_SITE_DONOR UTR_5_PRIME UTR_3_PRIME		
				do
					cat $INDEL_dir/$sampe.INDEL.tmp | grep "$indel" | cut -f2 | tr "," "\n" > $INDEL_dir/in.$sampe.$indel.tmp
					if [ ! -s $INDEL_dir/$sampe.$indel.tmp ]
					then
						echo "NOGENE"  >> $INDEL_dir/in.$sampe.$indel.tmp
					fi				
				done
				
				Rscript $script_path/summary.INDEL.r $report_dir/$sampe.gene.temp $INDEL_dir/in.$sampe.EXON_DELETED.tmp $INDEL_dir/in.$sampe.FRAME_SHIFT.tmp $INDEL_dir/in.$sampe.CODON_CHANGE.tmp $INDEL_dir/in.$sampe.UTR_5_DELETED.tmp $INDEL_dir/in.$sampe.UTR_3_DELETED.tmp $INDEL_dir/in.$sampe.CODON_INSERTION.tmp $INDEL_dir/in.$sampe.CODON_CHANGE_PLUS_CODON_INSERTION.tmp $INDEL_dir/in.$sampe.CODON_DELETION.tmp $INDEL_dir/in.$sampe.CODON_CHANGE_PLUS_CODON_DELETION.tmp $INDEL_dir/in.$sampe.SPLICE_SITE_ACCEPTOR.tmp $INDEL_dir/in.$sampe.SPLICE_SITE_DONOR.tmp $INDEL_dir/in.$sampe.UTR_5_PRIME.tmp $INDEL_dir/in.$sampe.UTR_3_PRIME.tmp $INDEL_dir/$sampe.EXON_DELETED.tmp $INDEL_dir/$sampe.FRAME_SHIFT.tmp $INDEL_dir/$sampe.CODON_CHANGE.tmp $INDEL_dir/$sampe.UTR_5_DELETED.tmp $INDEL_dir/$sampe.UTR_3_DELETED.tmp $INDEL_dir/$sampe.CODON_INSERTION.tmp $INDEL_dir/$sampe.CODON_CHANGE_PLUS_CODON_INSERTION.tmp $INDEL_dir/$sampe.CODON_DELETION.tmp $INDEL_dir/$sampe.CODON_CHANGE_PLUS_CODON_DELETION.tmp $INDEL_dir/$sampe.SPLICE_SITE_ACCEPTOR.tmp $INDEL_dir/$sampe.SPLICE_SITE_DONOR.tmp $INDEL_dir/$sampe.UTR_5_PRIME.tmp $INDEL_dir/$sampe.UTR_3_PRIME.tmp
				
				join $INDEL_dir/$sampe.EXON_DELETED.tmp $INDEL_dir/$sampe.FRAME_SHIFT.tmp > $INDEL_dir/$sampe.join1.txt
				join $INDEL_dir/$sampe.join1.txt $INDEL_dir/$sampe.CODON_CHANGE.tmp > $INDEL_dir/$sampe.join2.txt
				join $INDEL_dir/$sampe.join2.txt $INDEL_dir/$sampe.UTR_5_DELETED.tmp > $INDEL_dir/$sampe.join3.txt
				join $INDEL_dir/$sampe.join3.txt $INDEL_dir/$sampe.UTR_3_DELETED.tmp > $INDEL_dir/$sampe.join4.txt
				join $INDEL_dir/$sampe.join4.txt $INDEL_dir/$sampe.CODON_INSERTION.tmp > $INDEL_dir/$sampe.join5.txt
				join $INDEL_dir/$sampe.join5.txt $INDEL_dir/$sampe.CODON_CHANGE_PLUS_CODON_INSERTION.tmp > $INDEL_dir/$sampe.join6.txt
				join $INDEL_dir/$sampe.join6.txt $INDEL_dir/$sampe.CODON_DELETION.tmp > $INDEL_dir/$sampe.join7.txt
				join $INDEL_dir/$sampe.join7.txt $INDEL_dir/$sampe.CODON_CHANGE_PLUS_CODON_DELETION.tmp > $INDEL_dir/$sampe.join8.txt
				join $INDEL_dir/$sampe.join8.txt $INDEL_dir/$sampe.SPLICE_SITE_ACCEPTOR.tmp > $INDEL_dir/$sampe.join9.txt
				join $INDEL_dir/$sampe.join9.txt $INDEL_dir/$sampe.SPLICE_SITE_DONOR.tmp > $INDEL_dir/$sampe.join10.txt
				join $INDEL_dir/$sampe.join10.txt $INDEL_dir/$sampe.UTR_5_PRIME.tmp > $INDEL_dir/$sampe.join11.txt
				join $INDEL_dir/$sampe.join11.txt $INDEL_dir/$sampe.UTR_3_PRIME.tmp > $INDEL_dir/$sampe.join12.txt

				cat $INDEL_dir/$sampe.join12.txt | tr " " "\t" > $INDEL_dir/$sampe.join13.txt
				
				touch $INDEL_dir/$sampe.INDEL.summary
				echo -e "\tEXON_DELETED\tFRAME_SHIFT\tCODON_CHANGE\tUTR_5_DELETED\tUTR_3_DELETED\tCODON_INSERTION\tCODON_CHANGE_PLUS_CODON_INSERTION\tCODON_DELETION\tCODON_CHANGE_PLUS_CODON_DELETION\tSPLICE_SITE_ACCEPTOR\tSPLICE_SITE_DONOR\tUTR_5_PRIME\tUTR_3_PRIME" >> $INDEL_dir/$sampe.INDEL.summary
				cat $INDEL_dir/$sampe.join13.txt >> $INDEL_dir/$sampe.INDEL.summary
				Rscript $script_path/sum.cols.r $INDEL_dir/$sampe.INDEL.summary $INDEL_dir/$sampe.INDEL.sum
				
				rm $INDEL_dir/$sampe.*.tmp $INDEL_dir/in.$sampe.*.tmp $INDEL_dir/$sampe.join*.txt	
				
				
				cat $master_gene_file | awk '{print $4"\t"$1"\t"$2"\t"$3"\t"$5}' > $report_dir/$sampe.GeneList.forsummary.txt
				cat $master_entrez_file | awk '{print $2}' > $report_dir/$sampe.EntrezID.txt
				touch $report_dir/$sampe.Gene.Summary.txt
				echo -e "\t\t\t\t\t\t\t\t\t\t$sampe" >> $report_dir/$sampe.Gene.Summary.txt
				echo -e "\t\t\t\t\t\t\t\tSNV_BREAKDOWN\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tINDEL_BREAKDOWN" >> $report_dir/$sampe.Gene.Summary.txt
				echo -e "GENE\tCHROMOSOME\tSTART\tSTOP\tSTRAND\tENTREZ_GENE_ID\tTOTAL_SNVs\tTOTAL_INDELs\tSPLICE_SITE_ACCEPTOR\tSPLICE_SITE_DONOR\tSTART_LOST\tSTOP_GAINED\tSTOP_LOST\tRARE_AMINO_ACID\tNON_SYNONYMOUS_CODING\tSYNONYMOUS_START\tNON_SYNONYMOUS_START\tSTART_GAINED\tSYNONYMOUS_CODING\tSYNONYMOUS_STOP\tNON_SYNONYMOUS_STOP\tUTR_5_PRIME\tUTR_3_PRIME\tEXON_DELETED\tFRAME_SHIFT\tCODON_CHANGE\tUTR_5_DELETED\tUTR_3_DELETED\tCODON_INSERTION\tCODON_CHANGE_PLUS_CODON_INSERTION\tCODON_DELETION\tCODON_CHANGE_PLUS_CODON_DELETION\tSPLICE_SITE_ACCEPTOR\tSPLICE_SITE_DONOR\tUTR_5_PRIME\tUTR_3_PRIME" >> $report_dir/$sampe.Gene.Summary.txt
				sed -i '1d' $SNV_dir/$sampe.SNV.summary
				sed -i '1d' $INDEL_dir/$sampe.INDEL.summary
				cat $SNV_dir/$sampe.SNV.summary | cut -f2-16 > $SNV_dir/$sampe.SNV.tmp
				cat $INDEL_dir/$sampe.INDEL.summary | cut -f2-14 > $INDEL_dir/$sampe.INDEL.tmp
				
				paste $report_dir/$sampe.GeneList.forsummary.txt $report_dir/$sampe.EntrezID.txt $SNV_dir/$sampe.SNV.sum $INDEL_dir/$sampe.INDEL.sum $SNV_dir/$sampe.SNV.tmp $INDEL_dir/$sampe.INDEL.tmp >> $report_dir/$sampe.Gene.Summary.txt
				
				rm $SNV_dir/$sampe.SNV.tmp $INDEL_dir/$sampe.INDEL.tmp $report_dir/$sampe.GeneList.forsummary.txt $INDEL_dir/$sampe.INDEL.summary $INDEL_dir/$sampe.INDEL.sum $SNV_dir/$sampe.SNV.summary $SNV_dir/$sampe.SNV.sum $report_dir/$sampe.EntrezID.txt
				rm $report_dir/$sampe.gene.temp
			done				
		else
			### summarizing SNV files
			sample_info=$( cat $run_info | grep -w '^SAMPLE_INFO' | cut -d '=' -f2 )
			samples=$( cat $sample_info | grep -w "^$group" | cut -d '=' -f2 )
			let num_tumor=`echo $samples|tr " " "\n"|wc -l`-1
			tumor_list=`echo $samples | tr " " "\n" | tail -$num_tumor`
			if [ $somatic_calling == "NO" ]
			then
				tumor_list=$samples
			fi	
			for sample in $tumor_list    
			do
				if [ $somatic_calling == "YES" ]
				then
					sampe=$sampe
				else
					sampe=$sample
				fi	
				cat $master_gene_file | cut -f4 > $report_dir/$sampe.gene.temp
				if [ $somatic_calling == "NO" ]
				then
					file=$SNV_dir/$group.SNV.filtered.xls
				else
					file=$SNV_dir/TUMOR.$group.SNV.filtered.xls
				fi
				
				if [ ! -f $file ]
				then
					$script_path/errorlog.sh $file gene_summary.sh ERROR "not found"
					exit 1;
				fi	
				function=`awk -F '\t' '{ for(i=1;i<=NF;i++){ if ($i == "Effect") {print i} } }' $file`
				gene=`awk -F '\t' '{ for(i=1;i<=NF;i++){ if ($i == "geneList") {print i} } }' $file`
				sam=`awk -v sam=$sample -F '\t' '{ for(i=1;i<=NF;i++){ if ($i == sam) {print i} } }' $file`
				cat $file | awk 'NR>2' |  awk -v sam=$sam '$sam !~ /n\/a/'| cut -f "$function","$gene" > $SNV_dir/$sampe.SNV.tmp
				
				for snv in SPLICE_SITE_ACCEPTOR SPLICE_SITE_DONOR START_LOST STOP_GAINED STOP_LOST RARE_AMINO_ACID NON_SYNONYMOUS_CODING SYNONYMOUS_START NON_SYNONYMOUS_START START_GAINED SYNONYMOUS_CODING SYNONYMOUS_STOP NON_SYNONYMOUS_STOP UTR_5_PRIME UTR_3_PRIME
				do
					cat $SNV_dir/$sampe.SNV.tmp | grep "$snv" | cut -f2 | tr "," "\n" > $SNV_dir/in.$sampe.$snv.tmp
					if [ ! -s $SNV_dir/$sampe.$snv.tmp ]
					then
						echo "NOGENE"  >> $SNV_dir/in.$sampe.$snv.tmp
					fi	
				done    
			
				Rscript $script_path/summary.SNV.r $report_dir/$sampe.gene.temp $SNV_dir/in.$sampe.SPLICE_SITE_ACCEPTOR.tmp $SNV_dir/in.$sampe.SPLICE_SITE_DONOR.tmp $SNV_dir/in.$sampe.START_LOST.tmp $SNV_dir/in.$sampe.STOP_GAINED.tmp $SNV_dir/in.$sampe.STOP_LOST.tmp $SNV_dir/in.$sampe.RARE_AMINO_ACID.tmp $SNV_dir/in.$sampe.NON_SYNONYMOUS_CODING.tmp $SNV_dir/in.$sampe.SYNONYMOUS_START.tmp $SNV_dir/in.$sampe.NON_SYNONYMOUS_START.tmp $SNV_dir/in.$sampe.START_GAINED.tmp $SNV_dir/in.$sampe.SYNONYMOUS_CODING.tmp $SNV_dir/in.$sampe.SYNONYMOUS_STOP.tmp $SNV_dir/in.$sampe.NON_SYNONYMOUS_STOP.tmp $SNV_dir/in.$sampe.UTR_5_PRIME.tmp $SNV_dir/in.$sampe.UTR_3_PRIME.tmp $SNV_dir/$sampe.SPLICE_SITE_ACCEPTOR.tmp $SNV_dir/$sampe.SPLICE_SITE_DONOR.tmp $SNV_dir/$sampe.START_LOST.tmp $SNV_dir/$sampe.STOP_GAINED.tmp $SNV_dir/$sampe.STOP_LOST.tmp $SNV_dir/$sampe.RARE_AMINO_ACID.tmp $SNV_dir/$sampe.NON_SYNONYMOUS_CODING.tmp $SNV_dir/$sampe.SYNONYMOUS_START.tmp $SNV_dir/$sampe.NON_SYNONYMOUS_START.tmp $SNV_dir/$sampe.START_GAINED.tmp $SNV_dir/$sampe.SYNONYMOUS_CODING.tmp $SNV_dir/$sampe.SYNONYMOUS_STOP.tmp $SNV_dir/$sampe.NON_SYNONYMOUS_STOP.tmp $SNV_dir/$sampe.UTR_5_PRIME.tmp $SNV_dir/$sampe.UTR_3_PRIME.tmp
			
				join $SNV_dir/$sampe.SPLICE_SITE_ACCEPTOR.tmp $SNV_dir/$sampe.SPLICE_SITE_DONOR.tmp > $SNV_dir/$sampe.join1.txt
				join $SNV_dir/$sampe.join1.txt $SNV_dir/$sampe.START_LOST.tmp > $SNV_dir/$sampe.join2.txt
				join $SNV_dir/$sampe.join2.txt $SNV_dir/$sampe.STOP_GAINED.tmp > $SNV_dir/$sampe.join3.txt
				join $SNV_dir/$sampe.join3.txt $SNV_dir/$sampe.STOP_LOST.tmp > $SNV_dir/$sampe.join4.txt
				join $SNV_dir/$sampe.join4.txt $SNV_dir/$sampe.RARE_AMINO_ACID.tmp > $SNV_dir/$sampe.join5.txt
				join $SNV_dir/$sampe.join5.txt $SNV_dir/$sampe.NON_SYNONYMOUS_CODING.tmp > $SNV_dir/$sampe.join6.txt
				join $SNV_dir/$sampe.join6.txt $SNV_dir/$sampe.SYNONYMOUS_START.tmp > $SNV_dir/$sampe.join7.txt
				join $SNV_dir/$sampe.join7.txt $SNV_dir/$sampe.NON_SYNONYMOUS_START.tmp > $SNV_dir/$sampe.join8.txt
				join $SNV_dir/$sampe.join8.txt $SNV_dir/$sampe.START_GAINED.tmp > $SNV_dir/$sampe.join9.txt
				join $SNV_dir/$sampe.join9.txt $SNV_dir/$sampe.SYNONYMOUS_CODING.tmp > $SNV_dir/$sampe.join10.txt
				join $SNV_dir/$sampe.join10.txt $SNV_dir/$sampe.SYNONYMOUS_STOP.tmp > $SNV_dir/$sampe.join11.txt
				join $SNV_dir/$sampe.join11.txt $SNV_dir/$sampe.NON_SYNONYMOUS_STOP.tmp > $SNV_dir/$sampe.join12.txt
				join $SNV_dir/$sampe.join12.txt $SNV_dir/$sampe.UTR_5_PRIME.tmp > $SNV_dir/$sampe.join13.txt
				join $SNV_dir/$sampe.join13.txt $SNV_dir/$sampe.UTR_3_PRIME.tmp > $SNV_dir/$sampe.join14.txt
				cat $SNV_dir/$sampe.join14.txt | tr " " "\t" > $SNV_dir/$sampe.join15.txt
					
				touch $SNV_dir/$sampe.SNV.summary
				echo -e "\tSPLICE_SITE_ACCEPTOR\tSPLICE_SITE_DONOR\tSTART_LOST\tSTOP_GAINED\tSTOP_LOST\tRARE_AMINO_ACID\tNON_SYNONYMOUS_CODING\tSYNONYMOUS_START\tNON_SYNONYMOUS_START\tSTART_GAINED\tSYNONYMOUS_CODING\tSYNONYMOUS_STOP\tNON_SYNONYMOUS_STOP\tUTR_5_PRIME\tUTR_3_PRIME" >> $SNV_dir/$sampe.SNV.summary
				cat $SNV_dir/$sampe.join15.txt >> $SNV_dir/$sampe.SNV.summary
				Rscript $script_path/sum.cols.r $SNV_dir/$sampe.SNV.summary $SNV_dir/$sampe.SNV.sum
					
				rm $SNV_dir/$sampe.*.tmp $SNV_dir/in.$sampe.*.tmp $SNV_dir/$sampe.join*.txt

				#################################################################################################	
				### summarizing INDEL files           
				if [ $somatic_calling == "NO" ]
				then
					file=$INDEL_dir/$group.INDEL.filtered.xls
				else
					file=$INDEL_dir/TUMOR.$group.INDEL.filtered.xls
				fi	
				if [ ! -f $file ]
				then
					$script_path/errorlog.sh $file gene_summary.sh ERROR "not found"
					exit 1;
				fi	
				function=`awk -F '\t' '{ for(i=1;i<=NF;i++){ if ($i == "Effect") {print i} } }' $file`
				gene=`awk -F '\t' '{ for(i=1;i<=NF;i++){ if ($i == "geneList") {print i} } }' $file`
				sam=`awk -v sam=$sample -F '\t' '{ for(i=1;i<=NF;i++){ if ($i == sam) {print i} } }' $file`
				cat $file | awk 'NR>2' |  awk -v sam=$sam '$sam !~ /n\/a/' | cut -f "$function","$gene" > $INDEL_dir/$sampe.INDEL.tmp
				for indel in EXON_DELETED FRAME_SHIFT CODON_CHANGE UTR_5_DELETED UTR_3_DELETED CODON_INSERTION CODON_CHANGE_PLUS_CODON_INSERTION CODON_DELETION CODON_CHANGE_PLUS_CODON_DELETION SPLICE_SITE_ACCEPTOR SPLICE_SITE_DONOR UTR_5_PRIME UTR_3_PRIME		
				do
					cat $INDEL_dir/$sampe.INDEL.tmp | grep "$indel" | cut -f2 | tr "," "\n" > $INDEL_dir/in.$sampe.$indel.tmp
					if [ ! -s $INDEL_dir/$sampe.$indel.tmp ]
					then
						echo "NOGENE"  >> $INDEL_dir/in.$sampe.$indel.tmp
					fi				
				done
				
				Rscript $script_path/summary.INDEL.r $report_dir/$sampe.gene.temp $INDEL_dir/in.$sampe.EXON_DELETED.tmp $INDEL_dir/in.$sampe.FRAME_SHIFT.tmp $INDEL_dir/in.$sampe.CODON_CHANGE.tmp $INDEL_dir/in.$sampe.UTR_5_DELETED.tmp $INDEL_dir/in.$sampe.UTR_3_DELETED.tmp $INDEL_dir/in.$sampe.CODON_INSERTION.tmp $INDEL_dir/in.$sampe.CODON_CHANGE_PLUS_CODON_INSERTION.tmp $INDEL_dir/in.$sampe.CODON_DELETION.tmp $INDEL_dir/in.$sampe.CODON_CHANGE_PLUS_CODON_DELETION.tmp $INDEL_dir/in.$sampe.SPLICE_SITE_ACCEPTOR.tmp $INDEL_dir/in.$sampe.SPLICE_SITE_DONOR.tmp $INDEL_dir/in.$sampe.UTR_5_PRIME.tmp $INDEL_dir/in.$sampe.UTR_3_PRIME.tmp $INDEL_dir/$sampe.EXON_DELETED.tmp $INDEL_dir/$sampe.FRAME_SHIFT.tmp $INDEL_dir/$sampe.CODON_CHANGE.tmp $INDEL_dir/$sampe.UTR_5_DELETED.tmp $INDEL_dir/$sampe.UTR_3_DELETED.tmp $INDEL_dir/$sampe.CODON_INSERTION.tmp $INDEL_dir/$sampe.CODON_CHANGE_PLUS_CODON_INSERTION.tmp $INDEL_dir/$sampe.CODON_DELETION.tmp $INDEL_dir/$sampe.CODON_CHANGE_PLUS_CODON_DELETION.tmp $INDEL_dir/$sampe.SPLICE_SITE_ACCEPTOR.tmp $INDEL_dir/$sampe.SPLICE_SITE_DONOR.tmp $INDEL_dir/$sampe.UTR_5_PRIME.tmp $INDEL_dir/$sampe.UTR_3_PRIME.tmp
				
				join $INDEL_dir/$sampe.EXON_DELETED.tmp $INDEL_dir/$sampe.FRAME_SHIFT.tmp > $INDEL_dir/$sampe.join1.txt
				join $INDEL_dir/$sampe.join1.txt $INDEL_dir/$sampe.CODON_CHANGE.tmp > $INDEL_dir/$sampe.join2.txt
				join $INDEL_dir/$sampe.join2.txt $INDEL_dir/$sampe.UTR_5_DELETED.tmp > $INDEL_dir/$sampe.join3.txt
				join $INDEL_dir/$sampe.join3.txt $INDEL_dir/$sampe.UTR_3_DELETED.tmp > $INDEL_dir/$sampe.join4.txt
				join $INDEL_dir/$sampe.join4.txt $INDEL_dir/$sampe.CODON_INSERTION.tmp > $INDEL_dir/$sampe.join5.txt
				join $INDEL_dir/$sampe.join5.txt $INDEL_dir/$sampe.CODON_CHANGE_PLUS_CODON_INSERTION.tmp > $INDEL_dir/$sampe.join6.txt
				join $INDEL_dir/$sampe.join6.txt $INDEL_dir/$sampe.CODON_DELETION.tmp > $INDEL_dir/$sampe.join7.txt
				join $INDEL_dir/$sampe.join7.txt $INDEL_dir/$sampe.CODON_CHANGE_PLUS_CODON_DELETION.tmp > $INDEL_dir/$sampe.join8.txt
				join $INDEL_dir/$sampe.join8.txt $INDEL_dir/$sampe.SPLICE_SITE_ACCEPTOR.tmp > $INDEL_dir/$sampe.join9.txt
				join $INDEL_dir/$sampe.join9.txt $INDEL_dir/$sampe.SPLICE_SITE_DONOR.tmp > $INDEL_dir/$sampe.join10.txt
				join $INDEL_dir/$sampe.join10.txt $INDEL_dir/$sampe.UTR_5_PRIME.tmp > $INDEL_dir/$sampe.join11.txt
				join $INDEL_dir/$sampe.join11.txt $INDEL_dir/$sampe.UTR_3_PRIME.tmp > $INDEL_dir/$sampe.join12.txt

				cat $INDEL_dir/$sampe.join12.txt | tr " " "\t" > $INDEL_dir/$sampe.join13.txt
				
				touch $INDEL_dir/$sampe.INDEL.summary
				echo -e "\tEXON_DELETED\tFRAME_SHIFT\tCODON_CHANGE\tUTR_5_DELETED\tUTR_3_DELETED\tCODON_INSERTION\tCODON_CHANGE_PLUS_CODON_INSERTION\tCODON_DELETION\tCODON_CHANGE_PLUS_CODON_DELETION\tSPLICE_SITE_ACCEPTOR\tSPLICE_SITE_DONOR\tUTR_5_PRIME\tUTR_3_PRIME" >> $INDEL_dir/$sampe.INDEL.summary
				cat $INDEL_dir/$sampe.join13.txt >> $INDEL_dir/$sampe.INDEL.summary
				Rscript $script_path/sum.cols.r $INDEL_dir/$sampe.INDEL.summary $INDEL_dir/$sampe.INDEL.sum
				
				rm $INDEL_dir/$sampe.*.tmp $INDEL_dir/in.$sampe.*.tmp $INDEL_dir/$sampe.join*.txt	

				#################################################################################################	
				### summarizing CNV files
				if [ $somatic_calling == "NO" ]
				then
					sampe=$sample
				else
					sampe=$group.$sample
				fi	
				file=$CNV_dir/ANNOT/$sampe.CNV.annotated.txt
				gene=`awk -F '\t' '{ for(i=1;i<=NF;i++){ if ($i == "Gene") {print i} } }' $file`
				type=`awk -F '\t' '{ for(i=1;i<=NF;i++){ if ($i == "CNV_Type") {print i} } }' $file`
				cat $file | awk 'NR>2' | cut -f "$gene","$type" | tr "_" "\t" | cut -f1,3 | grep DUP > $CNV_dir/$sampe.DUP.tmp
				cat $file | awk 'NR>2' | cut -f "$gene","$type" | tr "_" "\t" | cut -f1,3 | grep DEL > $CNV_dir/$sampe.DEL.tmp

				if [ ! -s $CNV_dir/$sampe.DUP.tmp ]
				then
					echo "NOGENE" >> $CNV_dir/$sampe.DUP.tmp
				fi
				if [ ! -s $CNV_dir/$sampe.DEL.tmp ]
				then
					echo "NOGENE" >> $CNV_dir/$sampe.DEL.tmp
				fi

				Rscript $script_path/summary.CNV.r $report_dir/$sampe.gene.temp $CNV_dir/$sampe.DEL.tmp $CNV_dir/$sampe.DUP.tmp $CNV_dir/$sampe.DEL.txt $CNV_dir/$sampe.DUP.txt

				join $CNV_dir/$sampe.DEL.txt $CNV_dir/$sampe.DUP.txt > $CNV_dir/$sampe.join.txt
				cat $CNV_dir/$sampe.join.txt | tr " " "\t" > $CNV_dir/$sampe.join1.txt

				touch $CNV_dir/$sampe.CNV.summary
				echo -e "DEL\tDUP" >> $CNV_dir/$sampe.CNV.summary
				cat $CNV_dir/$sampe.join1.txt >> $CNV_dir/$sampe.CNV.summary
				Rscript $script_path/sum.cols.r $CNV_dir/$sampe.CNV.summary $CNV_dir/$sampe.CNV.sum

				rm $CNV_dir/$sampe.*.tmp $CNV_dir/$sampe.join*.txt $CNV_dir/$sampe.DEL.txt $CNV_dir/$sampe.DUP.txt
		#################################################################################################					### summarizing SV files
				file=$SV_dir/ANNOT/$sampe.SV.annotated.txt
				gene=`awk -F '\t' '{ for(i=1;i<=NF;i++){ if ($i == "GeneA_GeneB") {print i} } }' $file`
				type=`awk -F '\t' '{ for(i=1;i<=NF;i++){ if ($i == "SV_Type") {print i} } }' $file`
				for i in INV INS DEL ITX CTX
				do
					cat $file | awk 'NR>2' | cut -f "$type","$gene" | sed -e '/NA_/s//NOGENE_/g' -e '/_NA/s//_NOGENE/g' | tr "_" "\t" | cut -f2,3,4 | grep $i | tr " " "\t" > $SV_dir/$sampe.$i.tmp
				
					if [ ! -s $SV_dir/$sampe.$i.tmp ]
					then
						echo "NOGENE" >> $SV_dir/$sampe.$i.tmp
					fi
				done
				Rscript $script_path/summary.SV.r $report_dir/$sampe.gene.temp $SV_dir/$sampe.ITX.tmp $SV_dir/$sampe.INV.tmp $SV_dir/$sampe.DEL.tmp $SV_dir/$sampe.INS.tmp $SV_dir/$sampe.CTX.tmp $SV_dir/$sampe.ITX.txt $SV_dir/$sampe.INV.txt $SV_dir/$sampe.DEL.txt $SV_dir/$sampe.INS.txt $SV_dir/$sampe.CTX.txt
				
				join $SV_dir/$sampe.ITX.txt $SV_dir/$sampe.INV.txt > $SV_dir/$sampe.join1.txt
				join $SV_dir/$sampe.join1.txt $SV_dir/$sampe.DEL.txt > $SV_dir/$sampe.join2.txt
				join $SV_dir/$sampe.join2.txt $SV_dir/$sampe.INS.txt > $SV_dir/$sampe.join3.txt
				join $SV_dir/$sampe.join3.txt $SV_dir/$sampe.CTX.txt > $SV_dir/$sampe.join4.txt
				cat $SV_dir/$sampe.join4.txt | tr " " "\t" > $SV_dir/$sampe.join5.txt
				
				touch $SV_dir/$sampe.SV.summary
				echo -e "\tITX\tINV\tDEL\tINS\tCTX" >> $SV_dir/$sampe.SV.summary
				cat $SV_dir/$sampe.join5.txt >> $SV_dir/$sampe.SV.summary
				Rscript $script_path/sum.cols.r $SV_dir/$sampe.SV.summary $SV_dir/$sampe.SV.sum
				
				rm $SV_dir/$sampe.*.tmp $SV_dir/$sampe.join*.txt $SV_dir/$sampe.ITX.txt $SV_dir/$sampe.INV.txt $SV_dir/$sampe.DEL.txt $SV_dir/$sampe.INS.txt $SV_dir/$sampe.CTX.txt 
		#################################################################################################					### generating gene summary file
				cat $master_gene_file | awk '{print $4"\t"$1"\t"$2"\t"$3"\t"$5}' > $report_dir/$sampe.GeneList.forsummary.txt
				cat $master_entrez_file | awk '{print $2}' > $report_dir/$sampe.EntrezID.txt
				touch $report_dir/$sampe.Gene.Summary.txt
				echo -e "\t\t\t\t\t\t\t\t\t\t$sample" >> $report_dir/$sampe.Gene.Summary.txt
				echo -e "\t\t\t\t\t\t\t\t\t\tSNV_BREAKDOWN\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tINDEL_BREAKDOWN\t\t\t\t\t\t\t\t\t\t\t\t\tCNV_BREAKDOWN\t\tSV_BREAKDOWN\t\t\t\t" >> $report_dir/$sampe.Gene.Summary.txt
				echo -e "GENE\tCHROMOSOME\tSTART\tSTOP\tSTRAND\tENTREZ_GENE_ID\tTOTAL_SNVs\tTOTAL_INDELs\tTOTAL_CNVs\tTOTAL_SVs\tSPLICE_SITE_ACCEPTOR\tSPLICE_SITE_DONOR\tSTART_LOST\tSTOP_GAINED\tSTOP_LOST\tRARE_AMINO_ACID\tNON_SYNONYMOUS_CODING\tSYNONYMOUS_START\tNON_SYNONYMOUS_START\tSTART_GAINED\tSYNONYMOUS_CODING\tSYNONYMOUS_STOP\tNON_SYNONYMOUS_STOP\tUTR_5_PRIME\tUTR_3_PRIME\tEXON_DELETED\tFRAME_SHIFT\tCODON_CHANGE\tUTR_5_DELETED\tUTR_3_DELETED\tCODON_INSERTION\tCODON_CHANGE_PLUS_CODON_INSERTION\tCODON_DELETION\tCODON_CHANGE_PLUS_CODON_DELETION\tSPLICE_SITE_ACCEPTOR\tSPLICE_SITE_DONOR\tUTR_5_PRIME\tUTR_3_PRIME\tDELETION\tDUPLICATION\tITX\tINV\tDEL\tINS\tCTX" >> $report_dir/$sampe.Gene.Summary.txt
				sed -i '1d' $SNV_dir/$sampe.SNV.summary
				sed -i '1d' $INDEL_dir/$sampe.INDEL.summary
				sed -i '1d' $CNV_dir/$sampe.CNV.summary
				sed -i '1d' $SV_dir/$sampe.SV.summary
				cat $SNV_dir/$sampe.SNV.summary | cut -f 2-16 > $SNV_dir/$sampe.SNV.tmp
				cat $INDEL_dir/$sampe.INDEL.summary | cut -f 2-14 > $INDEL_dir/$sampe.INDEL.tmp
				cat $CNV_dir/$sampe.CNV.summary | cut -f 2,3 > $CNV_dir/$sampe.CNV.tmp
				cat $SV_dir/$sampe.SV.summary | cut -f 2-6 > $SV_dir/$sampe.SV.tmp
				
				paste $report_dir/$sampe.GeneList.forsummary.txt $report_dir/$sampe.EntrezID.txt $SNV_dir/$sampe.SNV.sum $INDEL_dir/$sampe.INDEL.sum $CNV_dir/$sampe.CNV.sum $SV_dir/$sampe.SV.sum $SNV_dir/$sampe.SNV.tmp $INDEL_dir/$sampe.INDEL.tmp $CNV_dir/$sampe.CNV.tmp $SV_dir/$sampe.SV.tmp >> $report_dir/$sampe.Gene.Summary.txt
				
				rm $SNV_dir/$sampe.SNV.tmp $INDEL_dir/$sampe.INDEL.tmp $CNV_dir/$sampe.CNV.tmp $SV_dir/$sampe.SV.tmp $report_dir/$sampe.GeneList.forsummary.txt $SV_dir/$sampe.SV.summary 
				rm $SV_dir/$sampe.SV.sum $CNV_dir/$sampe.CNV.summary $CNV_dir/$sampe.CNV.sum $INDEL_dir/$sampe.INDEL.summary $INDEL_dir/$sampe.INDEL.sum $SNV_dir/$sampe.SNV.summary $SNV_dir/$sampe.SNV.sum $report_dir/$sampe.EntrezID.txt
				rm $report_dir/$sampe.gene.temp
			done
		fi
	fi
	echo `date`
fi
	

