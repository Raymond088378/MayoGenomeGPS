#!/bin/bash
## baheti.saurabh@mayo.edu

shopt -s nocasematch

echo -e "\nPlease answer all the question in (yes/no) else directed {case insensitive}\n"

echo -e "Options: \n 1) Things to remember before running the workflow  \n 2)  Do the QC \n 3) Questions"

read -p "Enter your option(1 , 2 or 3) : " workflow
while [[ $workflow != "1" && $workflow != "2"  && $workflow != "3" ]]
do
	read -p "Enter your option(1 , 2 or 3) : " workflow
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
		echo -e "\n Available MOdules are :"
		echo -e "\n1. alignment"
		echo -e "\n2. external or mayo"
		echo -e "\n3. realignment or realign_mayo"
		echo -e "\n4. variant"
		echo -e "\n5. ontarget"
		echo -e "\n6. annotation"
		echo -e "\nMake sure you have correct capture kit and Ontarget file in the tool info file\n"
	else
		echo -e "\n Available Modules are :"
		echo -e "\n1. alignment"
		echo -e "\n2. external or mayo"
		echo -e "\n3. realignment or realign_mayo"
		echo -e "\n4. variant"
		echo -e "\n5. ontarget"
	fi
	
	read -p "Are these single samples or paired samples (single/paired) : " kind
	while [[ $kind != "paired" && $kind != "single" ]]
	do
		read -p "Are these single samples or paired samples (single/paired) : " kind
	done	
	if [[ $kind == "paired" ]]
	then
		echo -e "\n Available Modules are :"
		echo -e "\n2. external or mayo"
		echo -e "\n3. realignment or realign_mayo"
		echo -e "\n4. variant"
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


