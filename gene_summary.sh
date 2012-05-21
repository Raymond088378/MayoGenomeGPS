#!/bin/sh
	
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
    echo "\nUsage: </path/to/output directory> <sample name> </path/to/run_info.txt> </path/yo/Reports_per_Sample>";
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
    email=$( cat $run_info | grep -w '^EMAIL' | cut -d '=' -f2)
    type=$( cat $run_info | grep -w '^TYPE' | cut -d '=' -f2)
    script_path=$( cat $tool_info | grep -w '^WHOLEGENOME_PATH' | cut -d '=' -f2 )
    bedtools=$( cat $tool_info | grep -w '^BEDTOOLS' | cut -d '=' -f2 )
#	sample=$(cat $run_info | grep -w '^SAMPLENAMES' | cut -d '=' -f2 | tr ":" "\n" | head -n $SGE_TASK_ID | tail -n 1) 
#	group=$(cat $run_info | grep -w '^GROUPNAMES' | cut -d '=' -f2 | tr ":" "\n" | head -n $SGE_TASK_ID | tail -n 1) 
    master_gene_file=$( cat $tool_info | grep -w '^MASTER_GENE_FILE' | cut -d '=' -f2 )
    master_entrez_file=$( cat $tool_info | grep -w '^MASTER_ENTREZ_FILE' | cut -d '=' -f2 )
    email=$( cat $run_info | grep -w '^EMAIL' | cut -d '=' -f2)
    queue=$( cat $run_info | grep -w '^QUEUE' | cut -d '=' -f2)
    multi_sample=$( cat $run_info | grep -w '^MULTISAMPLE' | cut -d '=' -f2)

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
	    function=`awk -F '\t' '{ for(i=1;i<=NF;i++){ if ($i == "functionGVS") {print i} } }' $file`
	    gene=`awk -F '\t' '{ for(i=1;i<=NF;i++){ if ($i == "geneList") {print i} } }' $file`
		id=`awk -F '\t' '{ for(i=1;i<=NF;i++){ if ($i == "Entrez_id") {print i} } }' $file`
	    cat $file | awk 'NR>2' | cut -f "$function","$gene" > $SNV_dir/$sample.SNV.tmp	
		for snv in SPLICE_SITE_ACCEPTOR SPLICE_SITE_DONOR START_LOST STOP_GAINED STOP_LOST RARE_AMINO_ACID NON_SYNONYMOUS_CODING SYNONYMOUS_START NON_SYNONYMOUS_START START_GAINED SYNONYMOUS_CODING SYNONYMOUS_STOP NON_SYNONYMOUS_STOP UTR_5_PRIME UTR_3_PRIME
		do
			cat $SNV_dir/$sample.SNV.tmp | grep "$snv" | cut -f2 | tr "," "\n" > $SNV_dir/$sample.$snv.tmp
			if [ -s $SNV_dir/$sample.$snv.tmp ]
			then
				echo "NOGENE"  > $SNV_dir/$sample.$snv.tmp
			fi	
		done    
 
		Rscript $script_path/summary.SNV.r $report_dir/$sample.gene.temp $SNV_dir/$sample.SPLICE_SITE_ACCEPTOR.tmp $SNV_dir/$sample.SPLICE_SITE_DONOR.tmp $SNV_dir/$sample.START_LOST.tmp $SNV_dir/$sample.STOP_GAINED.tmp $SNV_dir/$sample.STOP_LOST.tmp $SNV_dir/$sample.RARE_AMINO_ACID.tmp $SNV_dir/$sample.NON_SYNONYMOUS_CODING.tmp $SNV_dir/$sample.SYNONYMOUS_START.tmp $SNV_dir/$sample.NON_SYNONYMOUS_START.txt $SNV_dir/$sample.START_GAINED.txt $SNV_dir/$sample.SYNONYMOUS_CODING.txt $SNV_dir/$sample.SYNONYMOUS_STOP.txt $SNV_dir/$sample.NON_SYNONYMOUS_STOP.txt $SNV_dir/$sample.UTR_5_PRIME.txt $SNV_dir/$sample.UTR_3_PRIME.txt
		
		join $SNV_dir/$sample.SPLICE_SITE_ACCEPTOR.txt $SNV_dir/$sample.SPLICE_SITE_DONOR.txt > $SNV_dir/$sample.join1.txt
		join $SNV_dir/$sample.join1.txt $SNV_dir/$sample.START_LOST.txt > $SNV_dir/$sample.join2.txt
		join $SNV_dir/$sample.join2.txt $SNV_dir/$sample.STOP_GAINED.txt > $SNV_dir/$sample.join3.txt
		join $SNV_dir/$sample.join3.txt $SNV_dir/$sample.STOP_LOST.txt > $SNV_dir/$sample.join4.txt
		join $SNV_dir/$sample.join4.txt $SNV_dir/$sample.RARE_AMINO_ACID.txt > $SNV_dir/$sample.join5.txt
		join $SNV_dir/$sample.join5.txt $SNV_dir/$sample.NON_SYNONYMOUS_CODING.txt > $SNV_dir/$sample.join6.txt
		join $SNV_dir/$sample.join6.txt $SNV_dir/$sample.SYNONYMOUS_START.txt > $SNV_dir/$sample.join7.txt
		join $SNV_dir/$sample.join7.txt $SNV_dir/$sample.NON_SYNONYMOUS_START.txt > $SNV_dir/$sample.join8.txt
		join $SNV_dir/$sample.join8.txt $SNV_dir/$sample.START_GAINED.txt > $SNV_dir/$sample.join9.txt
		join $SNV_dir/$sample.join9.txt $SNV_dir/$sample.SYNONYMOUS_CODING.txt > $SNV_dir/$sample.join10.txt
		join $SNV_dir/$sample.join10.txt $SNV_dir/$sample.SYNONYMOUS_STOP.txt > $SNV_dir/$sample.join11.txt
		join $SNV_dir/$sample.join11.txt $SNV_dir/$sample.NON_SYNONYMOUS_STOP.txt > $SNV_dir/$sample.join12.txt
		join $SNV_dir/$sample.join12.txt $SNV_dir/$sample.UTR_5_PRIME.txt > $SNV_dir/$sample.join13.txt
		join $SNV_dir/$sample.join13.txt $SNV_dir/$sample.UTR_3_PRIME.txt > $SNV_dir/$sample.join14.txt
		cat $SNV_dir/$sample.join14.txt | tr " " "\t" > $SNV_dir/$sample.join15.txt
			
		touch $SNV_dir/$sample.SNV.summary
		echo -e "SPLICE_SITE_ACCEPTOR\tSPLICE_SITE_DONOR\tSTART_LOST\tSTOP_GAINED\tSTOP_LOST\tRARE_AMINO_ACID\tNON_SYNONYMOUS_CODING\tSYNONYMOUS_START\tNON_SYNONYMOUS_START\tSTART_GAINED\tSYNONYMOUS_CODING\tSYNONYMOUS_STOP\tNON_SYNONYMOUS_STOP\tUTR_5_PRIME\tUTR_3_PRIME" >> $SNV_dir/$sample.SNV.summary
		cat $SNV_dir/$sample.join15.txt >> $SNV_dir/$sample.SNV.summary
			Rscript $script_path/sum.cols.r $SNV_dir/$sample.SNV.summary $SNV_dir/$sample.SNV.sum
			
			rm $SNV_dir/$sample.*.tmp $SNV_dir/$sample.join*.txt $SNV_dir/$sample.SPLICE_SITE_ACCEPTOR.txt $SNV_dir/$sample.SPLICE_SITE_DONOR.txt $SNV_dir/$sample.START_LOST.txt $SNV_dir/$sample.STOP_GAINED.txt $SNV_dir/$sample.STOP_LOST.txt $SNV_dir/$sample.RARE_AMINO_ACID.txt $SNV_dir/$sample.NON_SYNONYMOUS_CODING.txt $SNV_dir/$sample.SYNONYMOUS_START.txt $SNV_dir/$sample.NON_SYNONYMOUS_START.txt $SNV_dir/$sample.START_GAINED.txt $SNV_dir/$sample.SYNONYMOUS_CODING.txt $SNV_dir/$sample.SYNONYMOUS_STOP.txt $SNV_dir/$sample.UTR_5_PRIME.txt $SNV_dir/$sample.UTR_3_PRIME.txt 
	#################################################################################################	
			### summarizing INDEL files
            file=$INDEL_dir/$sample.INDEL.filtered.xls
            function=`awk -F '\t' '{ for(i=1;i<=NF;i++){ if ($i == "functionGVS") {print i} } }' $file`
			gene=`awk -F '\t' '{ for(i=1;i<=NF;i++){ if ($i == "geneList") {print i} } }' $file`
            cat $file | awk 'NR>2' | cut -f "$function","$gene" > $INDEL_dir/$sample.INDEL.tmp

            for indel in EXON_DELETED FRAME_SHIFT CODON_CHANGE UTR_5_DELETED UTR_3_DELETED CODON_INSERTION CODON_CHANGE_PLUS_CODON_INSERTION CODON_DELETION CODON_CHANGE_PLUS_CODON_DELETION SPLICE_SITE_ACCEPTOR SPLICE_SITE_DONOR UTR_5_PRIME UTR_3_PRIME		
            do
			cat $INDEL_dir/$sample.INDEL.tmp | grep "$indel" | cut -f2 | tr "," "\n" > $INDEL_dir/$sample.$indel.tmp
            done
			
			if [ -s $INDEL_dir/$sample.coding.tmp ]
			then
				echo "Not empty"
			else
				echo "NOGENE" >> $INDEL_dir/$sample.coding.tmp
			fi
			if [ -s $INDEL_dir/$sample.frameshift.tmp ]
			then
				echo "Not empty"
			else
				echo "NOGENE" >> $INDEL_dir/$sample.frameshift.tmp
			fi
			if [ -s $INDEL_dir/$sample.splice3.tmp ]
			then
				echo "Not empty"
			else
				echo "NOGENE" >> $INDEL_dir/$sample.splice3.tmp
			fi
			if [ -s $INDEL_dir/$sample.splice5.tmp ]
			then
				echo "Not empty"
			else
				echo "NOGENE" >> $INDEL_dir/$sample.splice5.tmp
			fi
			if [ -s $INDEL_dir/$sample.utr3.tmp ]
			then
				echo "Not empty"
			else
				echo "NOGENE" >> $INDEL_dir/$sample.utr3.tmp
			fi
			if [ -s $INDEL_dir/$sample.utr5.tmp ]
			then
				echo "Not empty"
			else
				echo "NOGENE" >> $INDEL_dir/$sample.utr5.tmp
			fi

			Rscript $script_path/summary.INDEL.r $report_dir/$sample.gene.temp $INDEL_dir/$sample.coding.tmp $INDEL_dir/$sample.frameshift.tmp $INDEL_dir/$sample.splice3.tmp $INDEL_dir/$sample.splice5.tmp $INDEL_dir/$sample.utr3.tmp $INDEL_dir/$sample.utr5.tmp $INDEL_dir/$sample.coding.txt $INDEL_dir/$sample.frameshift.txt $INDEL_dir/$sample.splice3.txt $INDEL_dir/$sample.splice5.txt $INDEL_dir/$sample.utr3.txt $INDEL_dir/$sample.utr5.txt
			
			join $INDEL_dir/$sample.coding.txt $INDEL_dir/$sample.frameshift.txt > $INDEL_dir/$sample.join1.txt
			join $INDEL_dir/$sample.join1.txt $INDEL_dir/$sample.splice3.txt > $INDEL_dir/$sample.join2.txt
			join $INDEL_dir/$sample.join2.txt $INDEL_dir/$sample.splice5.txt > $INDEL_dir/$sample.join3.txt
			join $INDEL_dir/$sample.join3.txt $INDEL_dir/$sample.utr3.txt > $INDEL_dir/$sample.join4.txt
			join $INDEL_dir/$sample.join4.txt $INDEL_dir/$sample.utr5.txt > $INDEL_dir/$sample.join5.txt
			cat $INDEL_dir/$sample.join5.txt | tr " " "\t" > $INDEL_dir/$sample.join6.txt
			
			touch $INDEL_dir/$sample.INDEL.summary
			echo -e "CODING\tFRAMESHIFT\tSPLICE-3\tSPLICE-5\tUTR-3\tUTR-5" >> $INDEL_dir/$sample.INDEL.summary
			cat $INDEL_dir/$sample.join6.txt >> $INDEL_dir/$sample.INDEL.summary
			Rscript $script_path/sum.cols.r $INDEL_dir/$sample.INDEL.summary $INDEL_dir/$sample.INDEL.sum
			
			rm $INDEL_dir/$sample.*.tmp $INDEL_dir/$sample.join*.txt $INDEL_dir/$sample.coding.txt $INDEL_dir/$sample.frameshift.txt $INDEL_dir/$sample.splice3.txt $INDEL_dir/$sample.splice5.txt $INDEL_dir/$sample.utr3.txt $INDEL_dir/$sample.utr5.txt	

#################################################################################################					### generating gene summary file
			cat $master_gene_file | awk '{print $4"\t"$1"\t"$2"\t"$3"\t"$5}' > $report_dir/$sample.GeneList.forsummary.txt
			cat $master_entrez_file | awk '{print $2}' > $report_dir/$sample.EntrezID.txt
			touch $report_dir/$sample.Gene.Summary.txt
			echo -e "\t\t\t\t\t\t\t\t\t\t$sample" >> $report_dir/$sample.Gene.Summary.txt
			echo -e "\t\t\t\t\t\t\t\t\t\tSNV_breakdown\t\t\t\t\t\t\t\tINDEL_breakdown\t\t\t\t\t\t" >> $report_dir/$sample.Gene.Summary.txt
			echo -e "Gene\tChromosome\tStart\tStop\tStrand\tEntrez_Gene_ID\tTotal_SNVs\tTotal_INDELs\tNONSENSE\tMISSENSE\tCODING-SYNONYMOUS\tCODING-NOTMOD3\tSPLICE-3\tSPLICE-5\tUTR-3\tUTR-5\tCODING\tFRAMESHIFT\tSPLICE-3\tSPLICE-5\tUTR-3\tUTR-5\t" >> $report_dir/$sample.Gene.Summary.txt
			sed -i '1d' $SNV_dir/$sample.SNV.summary
			sed -i '1d' $INDEL_dir/$sample.INDEL.summary
			cat $SNV_dir/$sample.SNV.summary | cut -f2,3,4,5,6,7,8,9 > $SNV_dir/$sample.SNV.tmp
			cat $INDEL_dir/$sample.INDEL.summary | cut -f2,3,4,5,6,7,8,9 > $INDEL_dir/$sample.INDEL.tmp
			
			paste $report_dir/$sample.GeneList.forsummary.txt $report_dir/$sample.EntrezID.txt $SNV_dir/$sample.SNV.sum $INDEL_dir/$sample.INDEL.sum $SNV_dir/$sample.SNV.tmp $INDEL_dir/$sample.INDEL.tmp >> $report_dir/$sample.Gene.Summary.txt
			
			rm $SNV_dir/$sample.SNV.tmp $INDEL_dir/$sample.INDEL.tmp $report_dir/$sample.GeneList.forsummary.txt $INDEL_dir/$sample.INDEL.summary $INDEL_dir/$sample.INDEL.sum $SNV_dir/$sample.SNV.summary $SNV_dir/$sample.SNV.sum $report_dir/$sample.EntrezID.txt
			rm $report_dir/$sample.gene.temp
			
#################################################################################################	

		else
			echo "Whole Genome Anlaysis"
			### summarizing SNV files
			file=$SNV_dir/$sample.SNV.cleaned_annot_filtered.xls
			function=`awk -F '\t' '{ for(i=1;i<=NF;i++){ if ($i == "functionGVS") {print i} } }' $file`
			gene=`awk -F '\t' '{ for(i=1;i<=NF;i++){ if ($i == "geneList") {print i} } }' $file`
			cat $file | cut -f "$function","$gene" > $SNV_dir/$sample.SNV.tmp
			cat $SNV_dir/$sample.SNV.tmp | grep nonsense | cut -f2 | tr "," "\n" > $SNV_dir/$sample.nonsense.tmp
			cat $SNV_dir/$sample.SNV.tmp | grep missense | cut -f2 | tr "," "\n" > $SNV_dir/$sample.missense.tmp
			cat $SNV_dir/$sample.SNV.tmp | grep coding-synonymous | cut -f2 | tr "," "\n" > $SNV_dir/$sample.codingsynonymous.tmp
			cat $SNV_dir/$sample.SNV.tmp | grep coding-notMod3 | cut -f2 | tr "," "\n" > $SNV_dir/$sample.codingnotMod3.tmp
			cat $SNV_dir/$sample.SNV.tmp | grep splice-3 | cut -f2 | tr "," "\n" > $SNV_dir/$sample.splice3.tmp
			cat $SNV_dir/$sample.SNV.tmp | grep splice-5 | cut -f2 | tr "," "\n" > $SNV_dir/$sample.splice5.tmp
			cat $SNV_dir/$sample.SNV.tmp | grep utr-3 | cut -f2 | tr "," "\n" > $SNV_dir/$sample.utr3.tmp
			cat $SNV_dir/$sample.SNV.tmp | grep utr-5 | cut -f2 | tr "," "\n" > $SNV_dir/$sample.utr5.tmp
			if [ -s $SNV_dir/$sample.nonsense.tmp ]
			then
				echo "Not empty"
			else
				echo "NOGENE" >> $SNV_dir/$sample.nonsense.tmp
			fi
			if [ -s $SNV_dir/$sample.missense.tmp ]
			then
				echo "Not empty"
			else
				echo "NOGENE" >> $SNV_dir/$sample.missense.tmp
			fi
			if [ -s $SNV_dir/$sample.codingsynonymous.tmp ]
			then
				echo "Not empty"
			else
				echo "NOGENE" >> $SNV_dir/$sample.codingsynonymous.tmp
			fi
			if [ -s $SNV_dir/$sample.codingnotMod3.tmp ]
			then
				echo "Not empty"
			else
				echo "NOGENE" >> $SNV_dir/$sample.codingnotMod3.tmp
			fi
			if [ -s $SNV_dir/$sample.splice3.tmp ]
			then
				echo "Not empty"
			else
				echo "NOGENE" >> $SNV_dir/$sample.splice3.tmp
			fi
			if [ -s $SNV_dir/$sample.splice5.tmp ]
			then
				echo "Not empty"
			else
				echo "NOGENE" >> $SNV_dir/$sample.splice5.tmp
			fi
			if [ -s $SNV_dir/$sample.utr3.tmp ]
			then
				echo "Not empty"
			else
				echo "NOGENE" >> $SNV_dir/$sample.utr3.tmp
			fi
			if [ -s $SNV_dir/$sample.utr5.tmp ]
			then
				echo "Not empty"
			else
				echo "NOGENE" >> $SNV_dir/$sample.utr5.tmp
			fi

			Rscript $script_path/summary.SNV.r $report_dir/$sample.gene.temp $SNV_dir/$sample.nonsense.tmp $SNV_dir/$sample.missense.tmp $SNV_dir/$sample.codingsynonymous.tmp $SNV_dir/$sample.codingnotMod3.tmp $SNV_dir/$sample.splice3.tmp $SNV_dir/$sample.splice5.tmp $SNV_dir/$sample.utr3.tmp $SNV_dir/$sample.utr5.tmp $SNV_dir/$sample.nonsense.txt $SNV_dir/$sample.missense.txt $SNV_dir/$sample.codingsynonymous.txt $SNV_dir/$sample.codingnotMod3.txt $SNV_dir/$sample.splice3.txt $SNV_dir/$sample.splice5.txt $SNV_dir/$sample.utr3.txt $SNV_dir/$sample.utr5.txt
			
			join $SNV_dir/$sample.nonsense.txt $SNV_dir/$sample.missense.txt > $SNV_dir/$sample.join1.txt
			join $SNV_dir/$sample.join1.txt $SNV_dir/$sample.codingsynonymous.txt > $SNV_dir/$sample.join2.txt
			join $SNV_dir/$sample.join2.txt $SNV_dir/$sample.codingnotMod3.txt > $SNV_dir/$sample.join3.txt
			join $SNV_dir/$sample.join3.txt $SNV_dir/$sample.splice3.txt > $SNV_dir/$sample.join4.txt
			join $SNV_dir/$sample.join4.txt $SNV_dir/$sample.splice5.txt > $SNV_dir/$sample.join5.txt
			join $SNV_dir/$sample.join5.txt $SNV_dir/$sample.utr3.txt > $SNV_dir/$sample.join6.txt
			join $SNV_dir/$sample.join6.txt $SNV_dir/$sample.utr5.txt > $SNV_dir/$sample.join7.txt
			cat $SNV_dir/$sample.join7.txt | tr " " "\t" > $SNV_dir/$sample.join8.txt
			
			touch $SNV_dir/$sample.SNV.summary
			echo -e "NONSENSE\tMISSENSE\tCODING-SYNONYMOUS\tCODING-NOTMOD3\tSPLICE-3\tSPLICE-5\tUTR-3\tUTR-5" >> $SNV_dir/$sample.SNV.summary
			cat $SNV_dir/$sample.join8.txt >> $SNV_dir/$sample.SNV.summary
			Rscript $script_path/sum.cols.r $SNV_dir/$sample.SNV.summary $SNV_dir/$sample.SNV.sum
			
			rm $SNV_dir/$sample.*.tmp $SNV_dir/$sample.join*.txt $SNV_dir/$sample.nonsense.txt $SNV_dir/$sample.missense.txt $SNV_dir/$sample.codingsynonymous.txt $SNV_dir/$sample.codingnotMod3.txt $SNV_dir/$sample.splice3.txt $SNV_dir/$sample.splice5.txt $SNV_dir/$sample.utr3.txt $SNV_dir/$sample.utr5.txt
	#################################################################################################	
			### summarizing INDEL files
			file=$INDEL_dir/$sample.INDEL.cleaned_annot_filtered.xls
			function=`awk -F '\t' '{ for(i=1;i<=NF;i++){ if ($i == "functionGVS") {print i} } }' $file`
			gene=`awk -F '\t' '{ for(i=1;i<=NF;i++){ if ($i == "geneList") {print i} } }' $file`
			cat $file |  cut -f "$function","$gene" > $INDEL_dir/$sample.INDEL.tmp
			cat $INDEL_dir/$sample.INDEL.tmp | grep coding | cut -f2 | tr "," "\n" > $INDEL_dir/$sample.coding.tmp
			cat $INDEL_dir/$sample.INDEL.tmp | grep frameshift | cut -f2 | tr "," "\n" > $INDEL_dir/$sample.frameshift.tmp
			cat $INDEL_dir/$sample.INDEL.tmp | grep splice-3 | cut -f2 | tr "," "\n" > $INDEL_dir/$sample.splice3.tmp
			cat $INDEL_dir/$sample.INDEL.tmp | grep splice-5 | cut -f2 | tr "," "\n" > $INDEL_dir/$sample.splice5.tmp
			cat $INDEL_dir/$sample.INDEL.tmp | grep utr-3 | cut -f2 | tr "," "\n" > $INDEL_dir/$sample.utr3.tmp
			cat $INDEL_dir/$sample.INDEL.tmp | grep utr-5 | cut -f2 | tr "," "\n" > $INDEL_dir/$sample.utr5.tmp
			if [ -s $INDEL_dir/$sample.coding.tmp ]
			then
				echo "Not empty"
			else
				echo "NOGENE" >> $INDEL_dir/$sample.coding.tmp
			fi
			if [ -s $INDEL_dir/$sample.frameshift.tmp ]
			then
				echo "Not empty"
			else
				echo "NOGENE" >> $INDEL_dir/$sample.frameshift.tmp
			fi
			if [ -s $INDEL_dir/$sample.splice3.tmp ]
			then
				echo "Not empty"
			else
				echo "NOGENE" >> $INDEL_dir/$sample.splice3.tmp
			fi
			if [ -s $INDEL_dir/$sample.splice5.tmp ]
			then
				echo "Not empty"
			else
				echo "NOGENE" >> $INDEL_dir/$sample.splice5.tmp
			fi
			if [ -s $INDEL_dir/$sample.utr3.tmp ]
			then
				echo "Not empty"
			else
				echo "NOGENE" >> $INDEL_dir/$sample.utr3.tmp
			fi
			if [ -s $INDEL_dir/$sample.utr5.tmp ]
			then
				echo "Not empty"
			else
				echo "NOGENE" >> $INDEL_dir/$sample.utr5.tmp
			fi

			Rscript $script_path/summary.INDEL.r $report_dir/$sample.gene.temp $INDEL_dir/$sample.coding.tmp $INDEL_dir/$sample.frameshift.tmp $INDEL_dir/$sample.splice3.tmp $INDEL_dir/$sample.splice5.tmp $INDEL_dir/$sample.utr3.tmp $INDEL_dir/$sample.utr5.tmp $INDEL_dir/$sample.coding.txt $INDEL_dir/$sample.frameshift.txt $INDEL_dir/$sample.splice3.txt $INDEL_dir/$sample.splice5.txt $INDEL_dir/$sample.utr3.txt $INDEL_dir/$sample.utr5.txt
			
			join $INDEL_dir/$sample.coding.txt $INDEL_dir/$sample.frameshift.txt > $INDEL_dir/$sample.join1.txt
			join $INDEL_dir/$sample.join1.txt $INDEL_dir/$sample.splice3.txt > $INDEL_dir/$sample.join2.txt
			join $INDEL_dir/$sample.join2.txt $INDEL_dir/$sample.splice5.txt > $INDEL_dir/$sample.join3.txt
			join $INDEL_dir/$sample.join3.txt $INDEL_dir/$sample.utr3.txt > $INDEL_dir/$sample.join4.txt
			join $INDEL_dir/$sample.join4.txt $INDEL_dir/$sample.utr5.txt > $INDEL_dir/$sample.join5.txt
			cat $INDEL_dir/$sample.join5.txt | tr " " "\t" > $INDEL_dir/$sample.join6.txt
			
			touch $INDEL_dir/$sample.INDEL.summary
			echo -e "CODING\tFRAMESHIFT\tSPLICE-3\tSPLICE-5\tUTR-3\tUTR-5" >> $INDEL_dir/$sample.INDEL.summary
			cat $INDEL_dir/$sample.join6.txt >> $INDEL_dir/$sample.INDEL.summary
			Rscript $script_path/sum.cols.r $INDEL_dir/$sample.INDEL.summary $INDEL_dir/$sample.INDEL.sum
			
			rm $INDEL_dir/$sample.*.tmp $INDEL_dir/$sample.join*.txt $INDEL_dir/$sample.coding.txt $INDEL_dir/$sample.frameshift.txt $INDEL_dir/$sample.splice3.txt $INDEL_dir/$sample.splice5.txt $INDEL_dir/$sample.utr3.txt $INDEL_dir/$sample.utr5.txt	

	#################################################################################################	
			### summarizing CNV files
			cat $CNV_dir/ANNOT/$sample.CNV.annotated.txt | tr "_" "\t" | cut -f4,10 | grep DUP > $CNV_dir/$sample.DUP.tmp
			cat $CNV_dir/ANNOT/$sample.CNV.annotated.txt | tr "_" "\t" | cut -f4,10 | grep DEL > $CNV_dir/$sample.DEL.tmp
			
			if [ -s $CNV_dir/$sample.DUP.tmp ]
			then
				echo "Not empty"
			else
				echo "NOGENE" >> $CNV_dir/$sample.DUP.tmp
			fi
			if [ -s $CNV_dir/$sample.DEL.tmp ]
			then
				echo "Not empty"
			else
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
			
			if [ -s $SV_dir/$sample.INV.tmp ]
			then
				echo "Not empty"
			else
				echo "NOGENE" >> $SV_dir/$sample.INV.tmp
			fi
			if [ -s $SV_dir/$sample.INS.tmp ]
			then
				echo "Not empty"
			else
				echo "NOGENE" >> $SV_dir/$sample.INS.tmp
			fi
			if [ -s $SV_dir/$sample.DEL.tmp ]
			then
				echo "Not empty"
			else
				echo "NOGENE" >> $SV_dir/$sample.DEL.tmp
			fi
			if [ -s $SV_dir/$sample.ITX.tmp ]
			then
				echo "Not empty"
			else
				echo "NOGENE" >> $SV_dir/$sample.ITX.tmp
			fi
			if [ -s $SV_dir/$sample.CTX.tmp ]
			then
				echo "Not empty"
			else
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
			echo -e "\t\t\t\t\t\t\t\t\t\tSNV_breakdown\t\t\t\t\t\t\t\tINDEL_breakdown\t\t\t\t\t\tCNV_breakdown\t\tSV_breakdown\t\t\t\t" >> $report_dir/$sample.Gene.Summary.txt
			echo -e "Gene\tChromosome\tStart\tStop\tStrand\tEntrez_Gene_ID\tTotal_SNVs\tTotal_INDELs\tTotal_CNVs\tTotal_SVs\tNONSENSE\tMISSENSE\tCODING-SYNONYMOUS\tCODING-NOTMOD3\tSPLICE-3\tSPLICE-5\tUTR-3\tUTR-5\tCODING\tFRAMESHIFT\tSPLICE-3\tSPLICE-5\tUTR-3\tUTR-5\tDELETION\tDUPLICATION\tITX\tINV\tDEL\tINS\tCTX" >> $report_dir/$sample.Gene.Summary.txt
			sed -i '1d' $SNV_dir/$sample.SNV.summary
			sed -i '1d' $INDEL_dir/$sample.INDEL.summary
			sed -i '1d' $CNV_dir/$sample.CNV.summary
			sed -i '1d' $SV_dir/$sample.SV.summary
			cat $SNV_dir/$sample.SNV.summary | cut -f2,3,4,5,6,7,8,9 > $SNV_dir/$sample.SNV.tmp
			cat $INDEL_dir/$sample.INDEL.summary | cut -f2,3,4,5,6,7,8,9 > $INDEL_dir/$sample.INDEL.tmp
			cat $CNV_dir/$sample.CNV.summary | cut -f2,3,4,5,6,7,8,9 > $CNV_dir/$sample.CNV.tmp
			cat $SV_dir/$sample.SV.summary | cut -f2,3,4,5,6,7,8,9 > $SV_dir/$sample.SV.tmp
			
			paste $report_dir/$sample.GeneList.forsummary.txt $report_dir/$sample.EntrezID.txt $SNV_dir/$sample.SNV.sum $INDEL_dir/$sample.INDEL.sum $CNV_dir/$sample.CNV.sum $SV_dir/$sample.SV.sum $SNV_dir/$sample.SNV.tmp $INDEL_dir/$sample.INDEL.tmp $CNV_dir/$sample.CNV.tmp $SV_dir/$sample.SV.tmp >> $report_dir/$sample.Gene.Summary.txt
			
			rm $SNV_dir/$sample.SNV.tmp $INDEL_dir/$sample.INDEL.tmp $CNV_dir/$sample.CNV.tmp $SV_dir/$sample.SV.tmp $report_dir/$sample.GeneList.forsummary.txt $SV_dir/$sample.SV.summary 
			rm $SV_dir/$sample.SV.sum $CNV_dir/$sample.CNV.summary $CNV_dir/$sample.CNV.sum $INDEL_dir/$sample.INDEL.summary $INDEL_dir/$sample.INDEL.sum $SNV_dir/$sample.SNV.summary $SNV_dir/$sample.SNV.sum $report_dir/$sample.EntrezID.txt
			rm $report_dir/$sample.gene.temp
		fi	
	else
		echo "Multi sample"
		### summarizing SNV files
		sample_info=$( cat $run_info | grep -w '^SAMPLE_INFO' | cut -d '=' -f2 )
		samples=$( cat $sample_info | grep -w "^$group" | cut -d '=' -f2 )
        let num_tumor=`echo $samples|tr " " "\n"|wc -l`-1
        tumor_list=`echo $samples | tr " " "\n" | tail -$num_tumor`
        for sample in $tumor_list    
        do
			cat $master_gene_file | cut -f4 > $report_dir/$group.$sample.gene.temp
			file=$SNV_dir/$group.$sample.SNV.cleaned_annot_filtered.xls
			function=`awk -F '\t' '{ for(i=1;i<=NF;i++){ if ($i == "functionGVS") {print i} } }' $file`
			gene=`awk -F '\t' '{ for(i=1;i<=NF;i++){ if ($i == "geneList") {print i} } }' $file`
			cat $file | cut -f "$function","$gene" > $SNV_dir/$group.$sample.SNV.tmp
			cat $SNV_dir/$group.$sample.SNV.tmp | grep nonsense | cut -f2 | tr "," "\n" > $SNV_dir/$group.$sample.nonsense.tmp
			cat $SNV_dir/$group.$sample.SNV.tmp | grep missense | cut -f2 | tr "," "\n" > $SNV_dir/$group.$sample.missense.tmp
			cat $SNV_dir/$group.$sample.SNV.tmp | grep coding-synonymous | cut -f2 | tr "," "\n" > $SNV_dir/$group.$sample.codingsynonymous.tmp
			cat $SNV_dir/$group.$sample.SNV.tmp | grep coding-notMod3 | cut -f2 | tr "," "\n" > $SNV_dir/$group.$sample.codingnotMod3.tmp
			cat $SNV_dir/$group.$sample.SNV.tmp | grep splice-3 | cut -f2 | tr "," "\n" > $SNV_dir/$group.$sample.splice3.tmp
			cat $SNV_dir/$group.$sample.SNV.tmp | grep splice-5 | cut -f2 | tr "," "\n" > $SNV_dir/$group.$sample.splice5.tmp
			cat $SNV_dir/$group.$sample.SNV.tmp | grep utr-3 | cut -f2 | tr "," "\n" > $SNV_dir/$group.$sample.utr3.tmp
			cat $SNV_dir/$group.$sample.SNV.tmp | grep utr-5 | cut -f2 | tr "," "\n" > $SNV_dir/$group.$sample.utr5.tmp
			if [ -s $SNV_dir/$group.$sample.nonsense.tmp ]
			then
				echo "Not empty"
			else
				echo "NOGENE" >> $SNV_dir/$group.$sample.nonsense.tmp
			fi
			if [ -s $SNV_dir/$group.$sample.missense.tmp ]
			then
				echo "Not empty"
			else
				echo "NOGENE" >> $SNV_dir/$group.$sample.missense.tmp
			fi
			if [ -s $SNV_dir/$group.$sample.codingsynonymous.tmp ]
			then
				echo "Not empty"
			else
				echo "NOGENE" >> $SNV_dir/$group.$sample.codingsynonymous.tmp
			fi
			if [ -s $SNV_dir/$group.$sample.codingnotMod3.tmp ]
			then
				echo "Not empty"
			else
				echo "NOGENE" >> $SNV_dir/$group.$sample.codingnotMod3.tmp
			fi
			if [ -s $SNV_dir/$group.$sample.splice3.tmp ]
			then
				echo "Not empty"
			else
				echo "NOGENE" >> $SNV_dir/$group.$sample.splice3.tmp
			fi
			if [ -s $SNV_dir/$group.$sample.splice5.tmp ]
			then
				echo "Not empty"
			else
				echo "NOGENE" >> $SNV_dir/$group.$sample.splice5.tmp
			fi
			if [ -s $SNV_dir/$group.$sample.utr3.tmp ]
			then
				echo "Not empty"
			else
				echo "NOGENE" >> $SNV_dir/$group.$sample.utr3.tmp
			fi
			if [ -s $SNV_dir/$group.$sample.utr5.tmp ]
			then
				echo "Not empty"
			else
				echo "NOGENE" >> $SNV_dir/$group.$sample.utr5.tmp
			fi

			Rscript $script_path/summary.SNV.r $report_dir/$group.$sample.gene.temp $SNV_dir/$group.$sample.nonsense.tmp $SNV_dir/$group.$sample.missense.tmp $SNV_dir/$group.$sample.codingsynonymous.tmp $SNV_dir/$group.$sample.codingnotMod3.tmp $SNV_dir/$group.$sample.splice3.tmp $SNV_dir/$group.$sample.splice5.tmp $SNV_dir/$group.$sample.utr3.tmp $SNV_dir/$group.$sample.utr5.tmp $SNV_dir/$group.$sample.nonsense.txt $SNV_dir/$group.$sample.missense.txt $SNV_dir/$group.$sample.codingsynonymous.txt $SNV_dir/$group.$sample.codingnotMod3.txt $SNV_dir/$group.$sample.splice3.txt $SNV_dir/$group.$sample.splice5.txt $SNV_dir/$group.$sample.utr3.txt $SNV_dir/$group.$sample.utr5.txt
			
			join $SNV_dir/$group.$sample.nonsense.txt $SNV_dir/$group.$sample.missense.txt > $SNV_dir/$group.$sample.join1.txt
			join $SNV_dir/$group.$sample.join1.txt $SNV_dir/$group.$sample.codingsynonymous.txt > $SNV_dir/$group.$sample.join2.txt
			join $SNV_dir/$group.$sample.join2.txt $SNV_dir/$group.$sample.codingnotMod3.txt > $SNV_dir/$group.$sample.join3.txt
			join $SNV_dir/$group.$sample.join3.txt $SNV_dir/$group.$sample.splice3.txt > $SNV_dir/$group.$sample.join4.txt
			join $SNV_dir/$group.$sample.join4.txt $SNV_dir/$group.$sample.splice5.txt > $SNV_dir/$group.$sample.join5.txt
			join $SNV_dir/$group.$sample.join5.txt $SNV_dir/$group.$sample.utr3.txt > $SNV_dir/$group.$sample.join6.txt
			join $SNV_dir/$group.$sample.join6.txt $SNV_dir/$group.$sample.utr5.txt > $SNV_dir/$group.$sample.join7.txt
			cat $SNV_dir/$group.$sample.join7.txt | tr " " "\t" > $SNV_dir/$group.$sample.join8.txt
			
			touch $SNV_dir/$group.$sample.SNV.summary
			echo -e "NONSENSE\tMISSENSE\tCODING-SYNONYMOUS\tCODING-NOTMOD3\tSPLICE-3\tSPLICE-5\tUTR-3\tUTR-5" >> $SNV_dir/$group.$sample.SNV.summary
			cat $SNV_dir/$group.$sample.join8.txt >> $SNV_dir/$group.$sample.SNV.summary
			Rscript $script_path/sum.cols.r $SNV_dir/$group.$sample.SNV.summary $SNV_dir/$group.$sample.SNV.sum
			
			rm $SNV_dir/$group.$sample.*.tmp $SNV_dir/$group.$sample.join*.txt $SNV_dir/$group.$sample.nonsense.txt $SNV_dir/$group.$sample.missense.txt $SNV_dir/$group.$sample.codingsynonymous.txt $SNV_dir/$group.$sample.codingnotMod3.txt $SNV_dir/$group.$sample.splice3.txt $SNV_dir/$group.$sample.splice5.txt $SNV_dir/$group.$sample.utr3.txt $SNV_dir/$group.$sample.utr5.txt
	#################################################################################################	
			### summarizing INDEL files
			file=$INDEL_dir/$group.$sample.INDEL.cleaned_annot_filtered.xls
			function=`awk -F '\t' '{ for(i=1;i<=NF;i++){ if ($i == "functionGVS") {print i} } }' $file`
			gene=`awk -F '\t' '{ for(i=1;i<=NF;i++){ if ($i == "geneList") {print i} } }' $file`
			cat $file |  cut -f "$function","$gene" > $INDEL_dir/$group.$sample.INDEL.tmp
			cat $INDEL_dir/$group.$sample.INDEL.tmp | grep coding | cut -f2 | tr "," "\n" > $INDEL_dir/$group.$sample.coding.tmp
			cat $INDEL_dir/$group.$sample.INDEL.tmp | grep frameshift | cut -f2 | tr "," "\n" > $INDEL_dir/$group.$sample.frameshift.tmp
			cat $INDEL_dir/$group.$sample.INDEL.tmp | grep splice-3 | cut -f2 | tr "," "\n" > $INDEL_dir/$group.$sample.splice3.tmp
			cat $INDEL_dir/$group.$sample.INDEL.tmp | grep splice-5 | cut -f2 | tr "," "\n" > $INDEL_dir/$group.$sample.splice5.tmp
			cat $INDEL_dir/$group.$sample.INDEL.tmp | grep utr-3 | cut -f2 | tr "," "\n" > $INDEL_dir/$group.$sample.utr3.tmp
			cat $INDEL_dir/$group.$sample.INDEL.tmp | grep utr-5 | cut -f2 | tr "," "\n" > $INDEL_dir/$group.$sample.utr5.tmp
			if [ -s $INDEL_dir/$group.$sample.coding.tmp ]
			then
				echo "Not empty"
			else
				echo "NOGENE" >> $INDEL_dir/$group.$sample.coding.tmp
			fi
			if [ -s $INDEL_dir/$group.$sample.frameshift.tmp ]
			then
				echo "Not empty"
			else
				echo "NOGENE" >> $INDEL_dir/$group.$sample.frameshift.tmp
			fi
			if [ -s $INDEL_dir/$group.$sample.splice3.tmp ]
			then
				echo "Not empty"
			else
				echo "NOGENE" >> $INDEL_dir/$group.$sample.splice3.tmp
			fi
			if [ -s $INDEL_dir/$group.$sample.splice5.tmp ]
			then
				echo "Not empty"
			else
				echo "NOGENE" >> $INDEL_dir/$group.$sample.splice5.tmp
			fi
			if [ -s $INDEL_dir/$group.$sample.utr3.tmp ]
			then
				echo "Not empty"
			else
				echo "NOGENE" >> $INDEL_dir/$group.$sample.utr3.tmp
			fi
			if [ -s $INDEL_dir/$group.$sample.utr5.tmp ]
			then
				echo "Not empty"
			else
				echo "NOGENE" >> $INDEL_dir/$group.$sample.utr5.tmp
			fi

			Rscript $script_path/summary.INDEL.r $report_dir/$group.$sample.gene.temp $INDEL_dir/$group.$sample.coding.tmp $INDEL_dir/$group.$sample.frameshift.tmp $INDEL_dir/$group.$sample.splice3.tmp $INDEL_dir/$group.$sample.splice5.tmp $INDEL_dir/$group.$sample.utr3.tmp $INDEL_dir/$group.$sample.utr5.tmp $INDEL_dir/$group.$sample.coding.txt $INDEL_dir/$group.$sample.frameshift.txt $INDEL_dir/$group.$sample.splice3.txt $INDEL_dir/$group.$sample.splice5.txt $INDEL_dir/$group.$sample.utr3.txt $INDEL_dir/$group.$sample.utr5.txt
			
			join $INDEL_dir/$group.$sample.coding.txt $INDEL_dir/$group.$sample.frameshift.txt > $INDEL_dir/$group.$sample.join1.txt
			join $INDEL_dir/$group.$sample.join1.txt $INDEL_dir/$group.$sample.splice3.txt > $INDEL_dir/$group.$sample.join2.txt
			join $INDEL_dir/$group.$sample.join2.txt $INDEL_dir/$group.$sample.splice5.txt > $INDEL_dir/$group.$sample.join3.txt
			join $INDEL_dir/$group.$sample.join3.txt $INDEL_dir/$group.$sample.utr3.txt > $INDEL_dir/$group.$sample.join4.txt
			join $INDEL_dir/$group.$sample.join4.txt $INDEL_dir/$group.$sample.utr5.txt > $INDEL_dir/$group.$sample.join5.txt
			cat $INDEL_dir/$group.$sample.join5.txt | tr " " "\t" > $INDEL_dir/$group.$sample.join6.txt
			
			touch $INDEL_dir/$group.$sample.INDEL.summary
			echo -e "CODING\tFRAMESHIFT\tSPLICE-3\tSPLICE-5\tUTR-3\tUTR-5" >> $INDEL_dir/$group.$sample.INDEL.summary
			cat $INDEL_dir/$group.$sample.join6.txt >> $INDEL_dir/$group.$sample.INDEL.summary
			Rscript $script_path/sum.cols.r $INDEL_dir/$group.$sample.INDEL.summary $INDEL_dir/$group.$sample.INDEL.sum
			
			rm $INDEL_dir/$group.$sample.*.tmp $INDEL_dir/$group.$sample.join*.txt $INDEL_dir/$group.$sample.coding.txt $INDEL_dir/$group.$sample.frameshift.txt $INDEL_dir/$group.$sample.splice3.txt $INDEL_dir/$group.$sample.splice5.txt $INDEL_dir/$group.$sample.utr3.txt $INDEL_dir/$group.$sample.utr5.txt	

	#################################################################################################	
			### summarizing CNV files
			cat $CNV_dir/ANNOT/$group.$sample.CNV.annotated.txt | tr "_" "\t" | cut -f4,10 | grep DUP > $CNV_dir/$group.$sample.DUP.tmp
			cat $CNV_dir/ANNOT/$group.$sample.CNV.annotated.txt | tr "_" "\t" | cut -f4,10 | grep DEL > $CNV_dir/$group.$sample.DEL.tmp
			
			if [ -s $CNV_dir/$group.$sample.DUP.tmp ]
			then
				echo "Not empty"
			else
				echo "NOGENE" >> $CNV_dir/$group.$sample.DUP.tmp
			fi
			if [ -s $CNV_dir/$group.$sample.DEL.tmp ]
			then
				echo "Not empty"
			else
				echo "NOGENE" >> $CNV_dir/$group.$sample.DEL.tmp
			fi

			Rscript $script_path/summary.CNV.r $report_dir/$group.$sample.gene.temp $CNV_dir/$group.$sample.DEL.tmp $CNV_dir/$group.$sample.DUP.tmp $CNV_dir/$group.$sample.DEL.txt $CNV_dir/$group.$sample.DUP.txt
			
			join $CNV_dir/$group.$sample.DEL.txt $CNV_dir/$group.$sample.DUP.txt > $CNV_dir/$group.$sample.join.txt
			cat $CNV_dir/$group.$sample.join.txt | tr " " "\t" > $CNV_dir/$group.$sample.join1.txt
			
			touch $CNV_dir/$group.$sample.CNV.summary
			echo -e "DEL\tDUP" >> $CNV_dir/$group.$sample.CNV.summary
			cat $CNV_dir/$group.$sample.join1.txt >> $CNV_dir/$group.$sample.CNV.summary
			Rscript $script_path/sum.cols.r $CNV_dir/$group.$sample.CNV.summary $CNV_dir/$group.$sample.CNV.sum
			
			rm $CNV_dir/$group.$sample.*.tmp $CNV_dir/$group.$sample.join*.txt $CNV_dir/$group.$sample.DEL.txt $CNV_dir/$group.$sample.DUP.txt
#################################################################################################					### summarizing SV files
			cat $SV_dir/ANNOT/$group.$sample.SV.annotated.txt | sed -e '/NA_/s//NOGENE_/g' -e '/_NA/s//_NOGENE/g' | grep INV | tr "_" "\t" | cut -f9,10,11 | tr " " "\t" > $SV_dir/$group.$sample.INV.tmp
			cat $SV_dir/ANNOT/$group.$sample.SV.annotated.txt | sed -e '/NA_/s//NOGENE_/g' -e '/_NA/s//_NOGENE/g' | grep INS | tr "_" "\t" | cut -f9,10,11 | tr " " "\t" > $SV_dir/$group.$sample.INS.tmp
			cat $SV_dir/ANNOT/$group.$sample.SV.annotated.txt | sed -e '/NA_/s//NOGENE_/g' -e '/_NA/s//_NOGENE/g' | grep DEL | tr "_" "\t" | cut -f9,10,11 | tr " " "\t" > $SV_dir/$group.$sample.DEL.tmp
			cat $SV_dir/ANNOT/$group.$sample.SV.annotated.txt | sed -e '/NA_/s//NOGENE_/g' -e '/_NA/s//_NOGENE/g' | grep ITX | tr "_" "\t" | cut -f9,10,11 | tr " " "\t" > $SV_dir/$group.$sample.ITX.tmp
			cat $SV_dir/ANNOT/$group.$sample.SV.annotated.txt | sed -e '/NA_/s//NOGENE_/g' -e '/_NA/s//_NOGENE/g' | grep CTX | tr "_" "\t" | cut -f9,10,11 | tr " " "\t" > $SV_dir/$group.$sample.CTX.tmp
			
			if [ -s $SV_dir/$group.$sample.INV.tmp ]
			then
				echo "Not empty"
			else
				echo "NOGENE" >> $SV_dir/$group.$sample.INV.tmp
			fi
			if [ -s $SV_dir/$group.$sample.INS.tmp ]
			then
				echo "Not empty"
			else
				echo "NOGENE" >> $SV_dir/$group.$sample.INS.tmp
			fi
			if [ -s $SV_dir/$group.$sample.DEL.tmp ]
			then
				echo "Not empty"
			else
				echo "NOGENE" >> $SV_dir/$group.$sample.DEL.tmp
			fi
			if [ -s $SV_dir/$group.$sample.ITX.tmp ]
			then
				echo "Not empty"
			else
				echo "NOGENE" >> $SV_dir/$group.$sample.ITX.tmp
			fi
			if [ -s $SV_dir/$group.$sample.CTX.tmp ]
			then
				echo "Not empty"
			else
				echo "NOGENE" >> $SV_dir/$group.$sample.CTX.tmp
			fi

			Rscript $script_path/summary.SV.r $report_dir/$group.$sample.gene.temp $SV_dir/$group.$sample.ITX.tmp $SV_dir/$group.$sample.INV.tmp $SV_dir/$group.$sample.DEL.tmp $SV_dir/$group.$sample.INS.tmp $SV_dir/$group.$sample.CTX.tmp $SV_dir/$group.$sample.ITX.txt $SV_dir/$group.$sample.INV.txt $SV_dir/$group.$sample.DEL.txt $SV_dir/$group.$sample.INS.txt $SV_dir/$group.$sample.CTX.txt
			
			join $SV_dir/$group.$sample.ITX.txt $SV_dir/$group.$sample.INV.txt > $SV_dir/$group.$sample.join1.txt
			join $SV_dir/$group.$sample.join1.txt $SV_dir/$group.$sample.DEL.txt > $SV_dir/$group.$sample.join2.txt
			join $SV_dir/$group.$sample.join2.txt $SV_dir/$group.$sample.INS.txt > $SV_dir/$group.$sample.join3.txt
			join $SV_dir/$group.$sample.join3.txt $SV_dir/$group.$sample.CTX.txt > $SV_dir/$group.$sample.join4.txt
			cat $SV_dir/$group.$sample.join4.txt | tr " " "\t" > $SV_dir/$group.$sample.join5.txt
			
			touch $SV_dir/$group.$sample.SV.summary
			echo -e "ITX\tINV\tDEL\tINS\tCTX" >> $SV_dir/$group.$sample.SV.summary
			cat $SV_dir/$group.$sample.join5.txt >> $SV_dir/$group.$sample.SV.summary
			Rscript $script_path/sum.cols.r $SV_dir/$group.$sample.SV.summary $SV_dir/$group.$sample.SV.sum
			
			rm $SV_dir/$group.$sample.*.tmp $SV_dir/$group.$sample.join*.txt $SV_dir/$group.$sample.ITX.txt $SV_dir/$group.$sample.INV.txt $SV_dir/$group.$sample.DEL.txt $SV_dir/$group.$sample.INS.txt $SV_dir/$group.$sample.CTX.txt 
#################################################################################################					### generating gene summary file
			cat $master_gene_file | awk '{print $4"\t"$1"\t"$2"\t"$3"\t"$5}' > $report_dir/$group.$sample.GeneList.forsummary.txt
			cat $master_entrez_file | awk '{print $2}' > $report_dir/$group.$sample.EntrezID.txt
			touch $report_dir/$group.$sample.Gene.Summary.txt
			echo -e "\t\t\t\t\t\t\t\t\t\t$sample" >> $report_dir/$group.$sample.Gene.Summary.txt
			echo -e "\t\t\t\t\t\t\t\t\t\tSNV_breakdown\t\t\t\t\t\t\t\tINDEL_breakdown\t\t\t\t\t\tCNV_breakdown\t\tSV_breakdown\t\t\t\t" >> $report_dir/$group.$sample.Gene.Summary.txt
			echo -e "Gene\tChromosome\tStart\tStop\tStrand\tEntrez_Gene_ID\tTotal_SNVs\tTotal_INDELs\tTotal_CNVs\tTotal_SVs\tNONSENSE\tMISSENSE\tCODING-SYNONYMOUS\tCODING-NOTMOD3\tSPLICE-3\tSPLICE-5\tUTR-3\tUTR-5\tCODING\tFRAMESHIFT\tSPLICE-3\tSPLICE-5\tUTR-3\tUTR-5\tDELETION\tDUPLICATION\tITX\tINV\tDEL\tINS\tCTX" >> $report_dir/$group.$sample.Gene.Summary.txt
			sed -i '1d' $SNV_dir/$group.$sample.SNV.summary
			sed -i '1d' $INDEL_dir/$group.$sample.INDEL.summary
			sed -i '1d' $CNV_dir/$group.$sample.CNV.summary
			sed -i '1d' $SV_dir/$group.$sample.SV.summary
			cat $SNV_dir/$group.$sample.SNV.summary | cut -f 2,3,4,5,6,7,8,9 > $SNV_dir/$group.$sample.SNV.tmp
			cat $INDEL_dir/$group.$sample.INDEL.summary | cut -f 2,3,4,5,6,7,8,9 > $INDEL_dir/$group.$sample.INDEL.tmp
			cat $CNV_dir/$group.$sample.CNV.summary | cut -f 2,3,4,5,6,7,8,9 > $CNV_dir/$group.$sample.CNV.tmp
			cat $SV_dir/$group.$sample.SV.summary | cut -f 2,3,4,5,6,7,8,9 > $SV_dir/$group.$sample.SV.tmp
			
			paste $report_dir/$group.$sample.GeneList.forsummary.txt $report_dir/$group.$sample.EntrezID.txt $SNV_dir/$group.$sample.SNV.sum $INDEL_dir/$group.$sample.INDEL.sum $CNV_dir/$group.$sample.CNV.sum $SV_dir/$group.$sample.SV.sum $SNV_dir/$group.$sample.SNV.tmp $INDEL_dir/$group.$sample.INDEL.tmp $CNV_dir/$group.$sample.CNV.tmp $SV_dir/$group.$sample.SV.tmp >> $report_dir/$group.$sample.Gene.Summary.txt
			
			rm $SNV_dir/$group.$sample.SNV.tmp $INDEL_dir/$group.$sample.INDEL.tmp $CNV_dir/$group.$sample.CNV.tmp $SV_dir/$group.$sample.SV.tmp $report_dir/$group.$sample.GeneList.forsummary.txt $SV_dir/$group.$sample.SV.summary 
			rm $SV_dir/$group.$sample.SV.sum $CNV_dir/$group.$sample.CNV.summary $CNV_dir/$group.$sample.CNV.sum $INDEL_dir/$group.$sample.INDEL.summary $INDEL_dir/$group.$sample.INDEL.sum $SNV_dir/$group.$sample.SNV.summary $SNV_dir/$group.$sample.SNV.sum $report_dir/$group.$sample.EntrezID.txt
			rm $report_dir/$group.$sample.gene.temp
		done
	fi
	echo `date`
fi
	

