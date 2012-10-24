#!/bin/bash
## baheti.saurabh@mayo.edu
## Oct 10th 2012
shopt -s nocasematch

echo -e "\nPlease answer all the question in (yes/no) else directed {case insensitive}\n"

echo -e "Options: \n 1) Things to remember before running the workflow \n 2) Information about the config files and help creating it. \n 3) Things to remember after running the workflow \n 4) Questions ?"

read -p "Enter your option(1 , 2, 3 or 4) : " workflow
while [[ $workflow != "1" && $workflow != "2"  && $workflow != "3"  && $workflow != "4" ]]
do
	read -p "Enter your option(1 , 2 , 3 or 4 ) : " workflow
done

### switch case starts

if [[ $workflow == "1" ]]
then
	echo -e "\nCopy the configuration files from the config folder of the current workflow version and make the changes according to the run and study"
	echo -e "\nCreate the four configuration files namely\n 1. Sample info file \n 2. Tool info file \n 3. Run info file \n 4. Memory info file \n"
	
	read -p "What kind of samples are you analyzing exome or whole_genome (exome/whole_genome) : " type
	while [[ $type != "exome" && $type != "whole_genome" ]]
	do
		read -p "What kind of samples are you analyzing exome or whole_genome (exome/whole_genome) : " type
	done	
	
	if [[ $type == "exome" ]]
	then
		echo -e "\n Available Modules are :\n1. alignment\n2. mayo or external\n3. realign-mayo or realignment\n4. variant\n5. ontarget\n6. annotation\nMake sure you have correct capture kit and Ontarget file in the tool info file\n"
	else
		echo -e "\n Available Modules are :\n1. alignment\n2. external or mayo\n3. realignment or realign-mayo\n4. variant\n5. ontarget"
	fi
	
	read -p "Are these single samples or paired samples (single/paired) : " kind
	while [[ $kind != "paired" && $kind != "single" ]]
	do
		read -p "Are these single samples or paired samples (single/paired) : " kind
	done	
	if [[ $kind == "paired" ]]
	then
		echo -e "\n Available Modules are :\n1. mayo or external\n2. realign-mayo or realignment\n3. variant"
	else 
		echo -e "\nAll the modules mentioned above are available"
	fi		

	echo -e "\nParameters you should be aware of"
	
	read -p "SNV_CALLER : " snv_caller
	if [[ $snv_caller == "GATK" || $snv_caller == "SNVMIX" || $snv_caller == "BEAUTY_EXOME" ]]
	then
		echo -e "\n SNV CALLER is set right"
	else
		echo -e "\n Incorrect SNV CALLER: avialable options are GATK , SNVMIX or BEAUTY_EXOME"
	fi
	
	if [[ $kind == "paired" ]]
	then
		read -p "SOMATIC_CALLER : " somatic_caller
		if [[ $somatic_caller == "SOMATICSNIPER" || $somatic_caller == "JOINTSNVMIX" || $somatic_caller == "MUTECT"  || $somatic_caller == "BEAUTY_EXOME" ]]
		then
			echo -e "\n SOMATIC CALLER is set right"
		else
			echo -e "\n Incorrect SOMATIC CALLER: available options are SOMATICSNIPER , JOINTSNVMIX , MUTECT or BEAUTY_EXOME"
		fi
	fi
	
	## FASTQC
	read -p "Did you looked at the FASTQC results : " fastqc
	while [[ $fastqc != "YES" && $fastqc != "NO" ]]
	do
		read -p "Did you looked at the FASTQC results : " fastqc
	done
	
	if [[ $fastqc == "YES" ]]
	then
		echo -e "\n ok, GOOD JOB !!"
	else
		echo -e "\nPlease review fastqc results:"
		echo "\\rcfcluster-cifs\data2\bsi\reports\<flowcell id>"
	fi

	read -p "Did you looked at the NGSPortal to validate the samples : " portal
	while [[ $portal != "YES" && $portal != "NO" ]]
	do
		read -p "Did you looked at the NGSPortal to validate the samples : " portal
	done
	
	if [[ $portal == "YES" ]]
	then
		echo -e "\n ok, GOOD JOB !!"
	else
		echo -e "\n Log into the Portal to validate the samples and indexes and lane information"
		echo " http://charlotte:8886/NGSPortal/ "
		echo -e "\nContact: \nDougherty, Gregory T. <Dougherty.Gregory@mayo.edu>\nHorton, Iain F. <Horton.Iain@mayo.edu>\nBockol, Matthew A. (Matt) <Bockol.Matthew@mayo.edu>"
	fi	
	
	read -p "Any/Further Questions : " ques
	while [[ $ques != "YES" && $ques != "NO" ]]
	do
		read -p "Any/Further Questions : " ques
	done	
	
	if [[ $ques == "NO" ]]
	then
		echo -e "\nok. you should be good to run the workflow"
	else
		echo -e "\nContact:\nRoss, Christian <ross.christian@mayo.edu>\nBaheti, Saurabh <baheti.saurabh@mayo.edu>"
	fi
elif [[ $workflow == "2" ]]
then
	## sample info file
	read -p "Do you want to know how to create sample info file ? : " sample_info
	while [[ $sample_info != "yes" && $sample_info != "no" ]]
	do
		read -p "Do you want to know how to create sample info file ? : " sample_info
	done
	if [[ $sample_info == "yes" ]]
	then
		echo -e "\nAll the example sample info files are located in the config folder of the workflow\n"
		read -p "what type of input files you have (FASTQ, BAM, VCF, TXT) : " file
		while [[ $file != "FASTQ" && $file != "BAM"  && $file != "VCF" && $file != "TXT" ]]
		do
			read -p "what type of input files you have (FASTQ, BAM, VCF, TXT)" file
		done
		if [[ $file == "FASTQ" ]]
		then
			echo -e "Important points: \n1. "FASTQ:" followed by Sample name follows '=' sign and then read1 and read2 are tab separated (specify the name of the FASTQ files after ‘=’. The name before ‘=’ is the short name for the sample) \n2. If you have multiple fastq for the sample then each group should tab seperated too. \n3. If you are running paired analysis then pair informaation should also come into this file. \n4. pair name follows '=' sign and then normal sample followed by as many tumor samples user has (all tab sepearted)\n"
		elif [[ $file == "BAM" ]]
		then
			echo -e "Important points: \n1. "BAM:" followed by sample name follows '=' sign and name of the bam file for that samples. \n2. If you have multiple bam file for a sample then all the bam files should be tab seperated on the same line. \n3.If you are running paired analysis then pair informaation should also come into this file. \n4. pair name follows '=' sign and then normal sample followed by as many tumor samples user has (all tab sepearted)\n"
		elif [[ $file == "VCF" ]]
		then
			read -p "what type of variants you have to annotate ? (BOTH,SNV, INDEL) : " variant
			while [[ $variant != "BOTH" && $variant != "SNV" && $variant != "INDEL" ]]
			do
				read -p "what type of variants you have to annotate ? (BOTH,SNV, INDEL) : " variant
			done
			echo -e "\nWe follow VCF 4.1 format, so make sure your VCF file follows the same\n"
			if [[ $varaint == "BOTH" ]]
			then
				echo -e "Important points: \n 1. User has both the snvs and indels in the vcf file. \n2. User can use same vcf file to point to get indels and snvs but need to specify the same file mulitple times. \n3.User should make sure that the vcf file has same sample name as they have in teh vcf header field. \n4. For snv calls "SNV:" followed by sample name folloes '=' sign and name of the vcf file. \n5. For indels calls "INDEL:" followed by sample name folloes '=' sign and name of the vcf file.\n"
			elif [[ $variant == "SNV" ]]
			then
				echo -e "Important points: \n 1. User has both the snvs and indels in the vcf file. \n2. User can use same vcf file to point to get indels and snvs but need to specify the same file mulitple times. \n3.User should make sure that the vcf file has same sample name as they have in teh vcf header field. \n4. For snv calls "SNV:" followed by sample name folloes '=' sign and name of the vcf file.\n"
			else
				echo -e "Important points: \n 1. User has both the snvs and indels in the vcf file. \n2. User can use same vcf file to point to get indels and snvs but need to specify the same file mulitple times. \n3.User should make sure that the vcf file has same sample name as they have in teh vcf header field. \n4. For indels calls "INDEL:" followed by sample name folloes '=' sign and name of the vcf file.\n"
			fi
		elif [[ $file == "TXT" ]]
		then
			read -p "what type of variants you have to annotate ? (BOTH,SNV, INDEL) : " variant
			while [[ $variant != "BOTH" && $variant != "SNV" && $variant != "INDEL" ]]
			do
				read -p "what type of variants you have to annotate ? (BOTH,SNV, INDEL) : " variant
			done
			echo -e "\nWe follow VCF 4.1 format coordinate system, so make sure your TXT file follows the same\n"
			if [[ $varaint == "BOTH" ]]
			then
				echo -e "Important points: \n 1. User has both the snvs and indels in the vcf file. \n2. For snv calls "SNV:" followed by sample name folloes '=' sign and name of the txt file. \n3. For indels calls "INDEL:" followed by sample name folloes '=' sign and name of the txt file. \n4. Both the files should have four columns namely Chromosome, position, Reference Allele, Alternate Allel. \n5. For Indels teh position is for the reference\n"
			elif [[ $variant == "SNV" ]]
			then
				echo -e "Important points: \n 1. User has both the snvs and indels in the vcf file. \n2. For snv calls "SNV:" followed by sample name folloes '=' sign and name of the txt file.\n3.the files should have four columns namely Chromosome, position, Reference Allele, Alternate Allel.\n"
				
			else
				echo -e "Important points: \n 1. User has both the snvs and indels in the vcf file. \n2. For indels calls "INDEL:" followed by sample name folloes '=' sign and name of the txt file.\n4. the files should have four columns namely Chromosome, position, Reference Allele, Alternate Allel. \n5. For Indels teh position is for the reference\n"
			fi
		fi		
 	else
		echo -e "\nOk, fine"
	fi		
	## runinfo file
	read -p "Do you want to know how to create run info file ? : " run_info
	while [[ $run_info != "yes" && $run_info != "no" ]]
	do
		read -p "Do you want to know how to create run info file ? : " run_info
	done
	if [[ $run_info == "YES" ]]
	then
		echo -e "\nAll the example run info files are located in the config folder of the workflow\n"
		echo -e "Important Paramters: \nANALYSIS\nFor Runs from upstairs lab you should specify realign-mayo if you are starting with bams  or mayo if you are starting with fastq's as we want to populate seconday dashboard data base \nFor Single sample exome/whole_genome: alignment,mayo,realign-mayo,external,realignment,ontarget,variant,annotation\nFor Paired Sample exome/whole_genome: alignment, mayo, external,realign-mayo,realignment,variant \nMULTISAMPLE\nYES if paired analsyis\nNO if single sample analsyis\nLABINDEXES\nIf its  multiplexing run then user should provide the index (index is always a 6 letter code for illumina data (don't include I in the index name)\nSNV_CALLER\nThere are three callers which user can use GATK, SNVMIX or BEAUTY_EXOME\nSOMATIC_CALLER\nIf you are doing Paired anlysis then user needs to specify the the caller from one of SOMATICSNIPER, JOINTSNVMIX, MUTECT or BEAUTY_EXOME\nALIGNER\nIf you are starting with FASTQ then user has an option to align the data uisng BWA or NOVOALIGN\nTOOL_INFO\nfull path to the tool info file\nSAMPLE_INFO\nFull path to the sample info file\nMEMORY_INFO\nFull path to the memory info file\n."
	else
		echo -e "\nOk, fine"
	fi	
	## tool info file
	read -p "Do you want to know how to create tool info file ? : " tool_info
	while [[ $tool_info != "yes" && $tool_info != "no" ]]
	do
		read -p "Do you want to know how to create tool info file ? : " tool_info
	done
	if [[ $tool_info == "YES" ]]
	then
		echo -e "\nAn example of tool info files is located in the config folder with all the default options for the workflow\n"
		echo -e "Important Paramters:"
		echo -e "\nTHREADS\t\t\tnumber of threads to use during alignment, variant calling etc. (4)
				\nREORDERSAM\t\tto make sure BAM file with a valid sequence dictionary. (NO)
				\nEMIT_ALL_SITES\t\tto call variant for all the positions in the given region. (NO)
				\nVARIANT_FILTER\t\tto filter the variant calls using VQSR. (YES)
				\nSOMATIC_VARIANT_FILTER\t\tto filter the somatic calls using VQSR. (YES)
				\nDEPTH_FILTER\t\tto filter the variant call using depth (0)
				\nTARGETTED\t\tto call the variants only on target region. (YES)
				\nMARKDUP\t\tto mark the duplicate reads in bam file	(YES)
				\nREMOVE_DUP\t\tto remove the duplidate reads from the bam file. (FALSE)
				\nREMOVE_ALIGNED_BAM\t\tto remove the aligned BAM files within the worflow execution. (YES)
				\nT_DEPTH_FILTER\t\tto filter the variant uisng total depth in multiple sample vcf (6)
				\nUPLOAD_TABLEBROWSER\t\tto upload the data to the Table Browser (YES)
				\nPLATFORM\t\tsequencing platform information	(illumina)
				\nCENTER\t\tsequencing center information (mayo)
				\nQUEUE\t\twhich queue to be used on RCF cluster (ngs-rand)
				\nSNVMIX2_params\t\tspeicfy commandline SNVMIX2 paramerts to be applied ()
				\nSNVMIX2_Filter\t\tspecify commandline snvmix filters (-p 0.8)
				\nUnifiedGenotyper_params\t\tspecify commandline UnifiedGenotyper parameters (-maxAlleles 5)
				\nSOMATIC_INDEL_params\t\tspecify commandline somatic indel parameters (--window_size 1000)
				\nSOMATIC_SNIPER_params\t\tspecify commandline somatic sniper paramerters (-q 20 -Q 20)
				\nMUTECT_params\t\tspecify commandline mutect paramters 
				\nBREAKDANCER_params\t\tspecify commandline breakdancer paramters (-c 5 -r 10)
				\nCREST_params\t\tspecify commandline crest paramters
				\nJSM_Filter\t\tspecify the filtering paramters for JSM (-prob 0.1)
				\nJOINTSNVMIX_params\t\tspecify commandline Joint SNVMIX paramters
				\nNOVO_params\t\tspecify the comamndline novolalign paramters (-g 60 -x 2 -i PE 425,80 -r Random --hdrhd off -v 120)
				\nBWA_params\t\tspecify commandline BWA paramters (-l 32 -t 4)
				\nVQSR_params_SNV\t\tspecify commandline VQSR for SNV paramters (-an QD -an HaplotypeScore -an MQRankSum -an ReadPosRankSum -an FS -an MQ -an DP --maxGaussians 4 --percentBadVariants 0.05)
				\nVQSR_params_INDEL\t\tspecify commandline VQSR for INDEL paramters (-an QD -an FS -an HaplotypeScore -an ReadPosRankSum --maxGaussians 4 --percentBadVariants 0.12 -std 10.0)
				\nPICARD_ReadGroup_params\t\tspecify commandline picard read group function paramters (PL=illumina CN=mayo LB=hg19 CREATE_INDEX=true)
				\nSNP_DISTANCE_INDEL\t\twindow size to look for a indel close to snp (10)
				\nREALIGN_params\t\tto skip the region if the reads are more than this (--maxReadsForRealignment 20000 --maxReadsInMemory 150000)	
				\nVCF_annotation_params\t\t features to annotate the vcf file (-A QualByDepth -A MappingQualityRankSumTest -A ReadPosRankSumTest -A HaplotypeScore -A DepthOfCoverage -A MappingQualityZero -A DepthPerAlleleBySample -A RMSMappingQuality -A FisherStrand -A ForwardReverseAlleleCounts )
				\nBLAT_params\t\tparameters used for blat querying of variant position (-w 50 -m 70 -t 90)
				\nCNVNATOR_BINSIZE\t\tbin size in CNVnator (1000)
				\nPCT_READS_SEGSEQ\t\tmin percent of reads for CNV (0.05)
				\nMINFOLD\t\tmin. fold (0.5)
				\nMAXFOLD\t\tmaximum fold (1.5)
				\nDISTGAP\t\tmin. distance from a gap (1000)
				\nBLAT_PORT\t\tport to use to initiate the blat server for CREST (50000)
				\nBLAT_SERVER\t\tserver for blat server (localhost)
				\nSTRUCT_DIST_GENE\t\tmax. distance to a gene(1000)
				\nSTRUCT_MIN_SUPPORT\t\tmin. reads to supoprt a SV uisng CREST (10)
				\nSTRUCT_MIN_IDENTITY\t\tmin. identity for CREST (0.9)
				\nSTRUCT_PCT_BLACKLIST\t\tpercent overlap of a cnv (1)
				\nMAX_FILE_HANDLE\t\tmax. intermediate/temporary files in a given instance of a script(100)
				\nMAX_READS_MEM_SORT\t\tmax. reads to keep in memory while sorting a bam file (2000000)
				\nTB_PORT\t\tport for Table Browser (8886)
				\nTB_HOST\t\thost for Table Browser (charlotte)
				\nJOB_LIMIT\t\tjob limit per user from RCF (3000)
				"
	else
		echo -e "\nOk, fine"
	fi	
	read -p "Do you want to know how to create memory info file ? : " memory_info
	while [[ $memory_info != "yes" && $memory_info != "no" ]]
	do
		read -p "Do you want to know how to create memory info file ? : " memory_info
	done
	if [[ $memory_info == "YES" ]]
	then
		echo "The memory requirements for each job and JVM are specified it and user should not touch this file unless the project needs it. If user need to make the changes in the file then you should have atleast extra pair of eyes to make sure nothing will break."
	else
		echo -e "\nOk, fine"
	fi	
	read -p "Do you have any other questions ? : " ques
	while [[ $ques != "yes" && $ques != "no" ]]
	do
		read -p "Do you have any other questions ? : " ques
	done
	if [[ $ques == "YES" ]]
	then
		echo -e "\nContact:\nRoss, Christian <ross.christian@mayo.edu>\nBaheti, Saurabh <baheti.saurabh@mayo.edu>"
	else
		echo -e "\nok. you should be good to run the workflow"
	fi	
elif [[ $workflow == "3" ]]
then 
	## Run the workflow
	read -p "Did you run the workflow : " run 
	while [[ $run != "YES" && $run != "NO" ]]
	do
		read -p "Did you run the workflow : " run
	done
	
	if [[ $run == "YES" ]]
	then
		echo -e "\n ok, GOOD JOB !!"
	else
		echo -e "\nCould you please Run the workflow and then come back to checklist again"
		exit 1;
	fi

	## QC
	read -p "Did you do the QC on the data : " qc
	while [[ $qc != "YES" && $qc != "NO" ]]
	do
		read -p "Did you do the QC on the data : " qc
	done  
	
	if [[ $qc == "YES" ]]
	then
		echo -e "\n ok, GOOD JOB !!"
	else
		echo -e "\nLook for QC information here : "
		echo "http://bioinformatics.mayo.edu/BMI/bin/view/Main/BioinformaticsCore/Analytics/PISupportLinksMainInternal"
	fi

	read -p "Did you run transfer and cleanup script : " clean
	while [[ $clean != "YES" && $clean != "NO" ]]
	do
		read -p "Did you run transfer and cleanup script : " clean
	done  
	
	if [[ $clean == "YES" ]]
	then
		echo -e "\n ok, GOOD JOB !!"
	else
		echo -e "\nLook for the clean up and transfer script in the scripts folder (transfer_clean.sh)"
	fi

	read -p "Did you Rsync the data from delivery to Isilon NL server : " rsync
	while [[ $rsync != "YES" && $rsync != "NO" ]]
	do
		read -p "Did you Rsync the data from delivery to Isilon NL server : " rsync
	done 
	
	if [[ $rsync == "YES" ]]
	then
		echo -e "\n ok, GOOD JOB !!"
	else 
		echo -e "\nGo To this webpage to the Rsync process : "
		echo "http://rcfcluster1.mayo.edu/cgi-cim/rsync-agtc.pl "
	fi


	read -p "Did you transfer the results to the windows share : " transfer
	while [[ $transfer != "YES" && $transfer != "NO" ]]
	do
		read -p "Did you transfer the results to the windows share : " transfer
	done 
	
	if [[ $transfer == "YES" ]]
	then
		echo -e "\n ok, GOOD JOB !!"
	else
		echo -e "\nYou need to transfer the results to windows share : "
		echo "\\ressrv08\BIC_Projects\PI_Support_Projects\<PI> "
	fi

	read -p "Do you have more questions or problems with workflow : " ques
	while [[ $ques != "YES" && $ques != "NO" ]]
	do
		read -p "Do you have more questions or problems with workflow : " ques
	done 
	
	if [[ $ques == "NO" ]]
	then
		echo -e "\nok. you are good to deliver the results"
	else
		echo -e "\nContact:\nRoss, Christian <ross.christian@mayo.edu>\nBaheti, Saurabh <baheti.saurabh@mayo.edu>"
	fi
else
	echo -e "\nContact:\nRoss, Christian <ross.christian@mayo.edu>\nBaheti, Saurabh <baheti.saurabh@mayo.edu>"
fi	


