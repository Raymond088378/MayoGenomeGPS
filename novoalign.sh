#!/bin/bash

####
# novoalign.sh
# 4/24/2013
# run novoalign on a fastq (zipped or unzipped fastq)
# baheti.saurabh@mayo.edu
####

### function
function check_variable()	{
	message=$1
	if [[ "$2" == "" ]] 
	then 
		echo "$message is not set correctly."
		exit 1
	fi		
}	

if [ $# != 5 ]
then
	echo -e "\nscript to run novoalign on a fastq \
	\nUsage: ./novoalign.sh </path/to/input fastq> <reads ':' seperated read1:read2> \
		<sample name> </path/to/outputbam> <tool info file>"
	exit 1;
fi
START=$(date +%s)
input=$1
reads=$2
sample=$3
outputbam=$4
tool_info=$5

if [ `echo $reads | tr ":" "\n" | wc -l` == 2 ]
then
	read1=`echo $reads| cut -f1 -d ':'`
	read2=`echo $reads| cut -f2 -d ':'`
	fastq="-f $input/$read1 $input/$read2"
else
	read1=$reads
	fastq="-f $input/$read1"
fi
			

### local variables
script_path=$( cat $tool_info | grep -w '^WORKFLOW_PATH' | cut -d '=' -f2 )
check_variable "$tool_info:WORKFLOW_PATH" $script_path
paramaters=$( cat $tool_info | grep -w '^NOVO_params' | cut -d '=' -f2)
check_variable "$tool_info:NOVO_params" $paramaters
novoalign=$( cat $tool_info | grep -w '^NOVOALIGN' | cut -d '=' -f2)
check_variable "$tool_info:NOVOALIGN" $novoalign
indexedgenome=$( cat $tool_info | grep -w '^NOVO_REF' | cut -d '=' -f2)
check_variable "$tool_info:NOVO_REF" $indexedgenome
center=$( cat $tool_info | grep -w '^CENTER' | cut -d '=' -f2 )
check_variable "$tool_info:CENTER" $center
platform=$( cat $tool_info | grep -w '^PLATFORM' | cut -d '=' -f2 )
check_variable "$tool_info:PLATFORM" $platform
samtools=$( cat $tool_info | grep -w '^SAMTOOLS' | cut -d '=' -f2)
check_variable "$tool_info:SAMTOOLS" $samtools

    
### checking the quality string
score=`perl $script_path/checkFastqQualityScores.pl $input/$read1 10000`
if [ $score -gt 65 ] 
then
	qual="-F ILMFQ"
else
	qual="-F STDFQ"
fi

### run novoalign 
$novoalign $paramaters -d $indexedgenome $qual $fastq \
            -o SAM "@RG\tID:$sample\tSM:$sample\tLB:$sample\tPL:$platform\tCN:$center" \
            	|  $samtools/samtools view -bS - > $outputbam
            	
### error checking
$samtools/samtools view -H $outputbam 1>$outputbam.header 2> $outputbam.novoalign.fix.log
if [[ `cat $outputbam.novoalign.fix.log | wc -l` -gt 0 || `cat $outputbam.header | wc -l` -le 0 ]]	
then
    $script_path/errorlog.sh $outputbam align_novo.sh ERROR "truncated or corrupt"
    exit 1;
else
    rm $outputbam.novoalign.fix.log
    echo $fastq | sed -e 's/-f//g' | xargs -t rm -Rf
fi	
rm $outputbam.header

END=$(date +%s)
DIFF=$(( $END - $START ))
echo "novoalign for $sample took $DIFF seconds"            	            	
            	
            
            
            