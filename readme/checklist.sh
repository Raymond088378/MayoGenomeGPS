#!/bin/bash
## baheti.saurabh@mayo.edu

shopt -s nocasematch

echo -e "\nPlease answer all the question in (yes/no) else directed {case insensitive}\n"

echo -e "Options: \n 1) Things to remember before running the workflow \n 2) Example config files \n 3) Do the QC \n 4) Questions ?"

read -p "Enter your option(1 , 2, 3 or 4) : " workflow
while [[ $workflow != "1" && $workflow != "2"  && $workflow != "3"  && $workflow != "4" ]]
do
	read -p "Enter your option(1 , 2 , 3 or 4 ) : " workflow
done

### switch case starts

if [[ $workflow == "1" ]]
then
	echo -e "\nCopy the configuration files from the config folder of the current workflow version and make the changes according to the run and study"
	echo -e "\nCreate the threee configuration files namely\n 1. Sample info file \n 2. Tool info file \n 3. Run info file"
	
	read -p "What kind of samples are you analyzing exome or whole_genome (exome/whole_genome) : " type
	while [[ $type != "exome" && $type != "whole_genome" ]]
	do
		read -p "What kind of samples are you analyzing exome or whole_genome (exome/whole_genome) : " type
	done	
	
	if [[ $type == "exome" ]]
	then
		echo -e "\n Available MOdules are :\n1. alignment\n2. external or mayo\n3. realignment or realign_mayo\n4. variant\n5. ontarget\n6. annotation\nMake sure you have correct capture kit and Ontarget file in the tool info file\n"
	else
		echo -e "\n Available Modules are :\n1. alignment\n2. external or mayo\n3. realignment or realign_mayo\n4. variant\n5. ontarget"
	fi
	
	read -p "Are these single samples or paired samples (single/paired) : " kind
	while [[ $kind != "paired" && $kind != "single" ]]
	do
		read -p "Are these single samples or paired samples (single/paired) : " kind
	done	
	if [[ $kind == "paired" ]]
	then
		echo -e "\n Available Modules are :\n2. external or mayo\n3. realignment or realign_mayo\n4. variant"
	else 
		echo -e "\nAll the modules mentioned above are available"
	fi		

	echo -e "\nParameters you should be aware of"
	
	read -p "SNV_CALLER : " snv_caller
	if [[ $snv_caller == "GATK" || $snv_caller == "SNVMIX" ]]
	then
		echo -e "\n SNV CALLER is set right"
	else
		echo -e "\n Incorrect SNV CALLER: avialable options are GATK or SNVMIX"
	fi
	
	if [[ $kind == "paired" ]]
	then
		read -p "SOMATIC_CALLER : " somatic_caller
		if [[ $somatic_caller == "SOMATICSNIPER" || $somatic_caller == "JOINTSNVMIX" || $soamtic_caller == "MUTECT" ]]
		then
			echo -e "\n SOMATIC CALLER is set right"
		else
			echo -e "\n Incorrect SOMATIC CALLER: available options are SOMATICSNIPER , JOINTSNVMIX or MUTECT"
		fi
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
		echo -e "\nContact: \n Baheti, Saurabh <baheti.saurabh@mayo.edu> \n Davila, Jaime I., Ph.D. <Davila.Jaime@mayo.edu>"
	fi
elif [[ $workflow == "2" ]]
then
	## sample info file
	read -p "Do you want to know how to create sample info file ? : " sample_info
	while [[ $sample_info != "yes" && $sample_info != "no" ]]
	do
		read -p "Do you want to know how to create sample info file ? : " sample_info
	done
	if [[ $sample_info == "YES" ]]
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
		echo -e "Important Paramters: \nANALYSIS\nFor Runs from upstairs lab you should specify realign-mayo if you are starting with bams  or mayo if you are starting with fastq's as we want to populate seconday dashboard data base \nFor Single sample exome/whole_genome: alignment,mayo,realign-mayo,external,realignment,ontarget,variant,annotation\nFor Paired Sample exome/whole_genome: alignment, mayo, external,realign-mayo,realignment,variant \nMULTISAMPLE\nYES if paired analsyis\nNO if single sample analsyis\nLABINDEXES\nIf its  multiplexing run then user should provide the index (index is always a 6 letter code for illumina data (don't include I in the index name)\nSNV_CALLER\nThere are two callers which user can use GATK or SNVMIX\nSOMATIC_CALLER\nIf you are doing Paired anlysis then user needs to specify the the caller from one of SOMATICSNIPER, JOINTSNVMIX or MUTECT\nUPLOAD_TABLEBROWSER\nthis should always be YES for all the mayo runs\nALIGNER\nIf you are starting with FASTQ then user has an option to align the data uisng BWA or NOVOALIGN.\n"
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
		echo -e "\nSOMATIC_THRESHOLD				used for joint snvmix to have a cutoff at min 0.1 as somatic probability (0.1)
				\nT_DEPTH_FILTER					min number of supported reads for multi sample calling (6)
				\nSOMATIC_QUALITY					min somatic qaulity required for somatic snipper to call variants (40)
				\nMAPPING_QUALITY					min mapping quality required to call variants (20)
				\nBASE_QUALITY						min base quality required to call the variants (20)
				\nEMIT_ALL_SITES					Boolean (YES/NO) Yes if user needs to call variant at all the postiions, NO otherwise (NO)
				\nVARIANT_FILTER					If user wants to filter the variant using our recommentations (YES)
				\nDEPTH_FILTER						If user wants to filter the variant using a depth filtering (0)
				\nTARGETTED							If user wants to call variant only in the target region, normal practice is to call teh variant on whole gonome and then do the subset (NO)
				\nREMOVE_DUP						If user want to flag or remove the duplicate reads (FALSE)
				\nPROB_FILTER						filter the variant using probability threshold for snvmix (0.8)
				\nSOMATIC_INDEL_FILTER_EXPRESSION	somatic indel calling filter expression (NA)
				\nMIN_SCIP_LEN						minmum length of soft clipped reads required to call SV using CREST(12)
				\nMIN_SCIP_READS					min. number of soft clipped reads required to call SV uising CREST (14)
				\nMAX_ALT_ALLELES					Max alternate allele to report during genotyping (5)
				\nSNP_DISTANCE_INDEL				window size to look for a indel close to snp (10)
				\nWINDOW_BLAT						window size to blat a position to teh refernce genome (50)
				\n
				"
	fi	
	read -p "Do you have any other questions ? : " ques
	while [[ $ques != "yes" && $ques != "no" ]]
	do
		read -p "Do you have any other questions ? : " ques
	done
	if [[ $ques == "YES" ]]
	then
		echo -e "\nContact: \n Baheti, Saurabh <baheti.saurabh@mayo.edu> \n Davila, Jaime I., Ph.D. <Davila.Jaime@mayo.edu>"
	else
		echo -e "\nok. you should be good to run the workflow"
	fi	
elif [[ $workflow == "3" ]]
then 
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
		echo -e "\nPlease review fastqc results here"
		echo "\\rcfcluster-cifs\data2\bsi\reports\<flowcell id>"
		exit 1;
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
		exit 1;
	fi	

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
		exit 1;
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
		exit 1;
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
		echo -e "\nGo To this webpage to the R sync process : "
		echo "http://rcfcluster1.mayo.edu/cgi-cim/rsync-agtc.pl "
		exit 1;
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
		exit 1;
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
		echo -e "Contact: \n Baheti, Saurabh <baheti.saurabh@mayo.edu> \n Davila, Jaime I., Ph.D. <Davila.Jaime@mayo.edu>"
	fi
else
	echo -e "\nContact: \n Baheti, Saurabh <baheti.saurabh@mayo.edu> \n Davila, Jaime I., Ph.D. <Davila.Jaime@mayo.edu>"
fi	


