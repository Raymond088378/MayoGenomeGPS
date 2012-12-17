#! /bin/bash
### Written by: Raymond Moore
### Purpose: Stand-Alone IM Module for Pharmicogenetic annotation
### Input: Tab-Delinated file (likely Resulting from Treat Workflow)
### Output: Appended Phamicogenetic Annotation (tabbed) file
  
set -x
if [ "$#" -lt "2" ]
then	
	echo "Usage: <full path to treat file> <run info>";
else
	###################################
	####  SET UP VARIABLES TO RUN  ####
	###################################
	echo `date`
	INPUT=$1
	run_info=$2
	tool_info=$( cat $run_info | grep -w '^TOOL_INFO' | cut -d '=' -f2)
	java=$( cat $tool_info | grep -w '^JAVA' | cut -d '=' -f2 )
	DIR="$( cat $tool_info | grep -w '^WORKFLOW_PATH' | cut -d '=' -f2 )/beauty_annot_module"
	NUM=$( cat $tool_info | grep -w '^BEAUTY_MOD_MULTIPLIER' | cut -d '=' -f2 )
	DB=$( cat $tool_info | grep -w '^ANNOTATION_MODULE_DATA' | cut -d '=' -f2 )
	feature_selection=$( cat $tool_info | grep -w '^ANNOTATION_FEATURES' | cut -d '=' -f2 )
	### Depricated	
	#QUEUE=$( cat $tool_info | grep -w '^QUEUE' | sed -e '/QUEUE=/s///g')
	#EMAIL="here@mayo.edu"
	###
	output=$( cat $run_info | grep -w '^BASE_OUTPUT_DIR' | cut -d '=' -f2)
	PI=$( cat $run_info | grep -w '^PI' | cut -d '=' -f2)
	tool=$( cat $run_info | grep -w '^TYPE' | cut -d '=' -f2|tr "[A-Z]" "[a-z]")
	run_num=$( cat $run_info | grep -w '^OUTPUT_FOLDER' | cut -d '=' -f2)
	output_dir=$output/$PI/$tool/$run_num
	perl="/usr/local/biotools/perl/5.14.2/bin/perl"

	email=`finger $USER | grep Mail | awk '{print $NF}'`
	queue=$( cat $tool_info | grep -w '^QUEUE' | sed -e '/QUEUE=/s///g')
	ARGS="#\$ -V\n#\$ -cwd\n#\$ -o $output_dir/logs\n#\$ -e $output_dir/logs\n#\$ -q $queue\n#\$ -m a\n#\$ -M $email\n#\$ -l h_stack=10M"

	####
	# Check for input file
	####
	if [ ! -f $INPUT ];
	then
		echo "Treat Input File not found!"
		exit 1;
	fi
	BASENOM=${INPUT##*/} ### filename + extension (no path)
	INEXT=${BASENOM##*.} ### input extension (xls or tsv or txt)
	PREPATH=${INPUT%/*} ### provided path to input directory
	RAWNAME=${BASENOM%.*} ### filename no extension
	####
	# Check for Full Path
	####
	if [ ! -d "$PREPATH" ]
	then
		echo "Must Provide Full Path To File!"
		exit 1;
	fi
	####
	# Check for Feature Selections file/xls
	####
	default=$PREPATH"/"$RAWNAME"_FeatureSelection.xls"
 	FEATURE=${feature_selection:-$default}
	if [ ! -f $FEATURE ];
	then
		echo "Feature Selection File not found!"
		exit 1;
	fi
	
	TMPDIR="$PREPATH/tmp$RANDOM";
	mkdir -p $TMPDIR
	mkdir -p "$TMPDIR/scpts"
	LOGS="$TMPDIR/Logs";
	mkdir -p $LOGS

	####
	# Calculate number of lines per file, based on how many nodes requirested
	####
	LINES=`cat $INPUT | wc -l`;
	(( LINES += 1 ))
	SPLIT=$(($LINES/$NUM))
	
	##
	# Split File into sections to Array Iterate over.
	##
	SPLITTING=$(perl $DIR/splitTREAT.pl $INPUT $FEATURE $TMPDIR $SPLIT  2>&1)
	if [  "$SPLITTING" != "" ]
	then
		echo -e "ERROR:\n$SPLITTING"
		#rm -R $TMPDIR
		exit 1;
	fi
		
	##########################################
	####  PARALLEL ANNOTATION COLLECTION  ####
	##########################################
	
	######
	## Function to Set up & Run SGE.
	## SGE_JOB_ID=$(sendtoSGE <compile> <script> <options + delin> <vmax> <job hold>);
	######
	sendtoSGE(){
		### Expect 1:execution 2:script 3:options(-d"|") 4:vMem 5:ArrayJob? 6:hold(optional)
		### Output SGE job Id
		JBHLD="";
		if [  $# -eq 6 ]
		then
			JBHLD="\n#$ -hold_jid $6";
		fi
		if [ $5 == "1" ]
		then
			ARRAY="#\$ -t 1-$NUM:1\n"
		else
			ARRAY=""
		fi
		if [  $1 == "java" ]
		then
			EXE="$java/java -jar -Xms200m -Xmx1g"
		else
			EXE=$1
		fi
		scpt=${2%.*}
		TMPNAME="$TMPDIR/scpts/${scpt}_${RANDOM}.sh"
		SCRIPT="#! /bin/bash\n $ARGS\n${ARRAY}\n#\$ -l h_vmem=$4\n#\$ $JBHLD"
		echo -e $SCRIPT >> $TMPNAME
		OPTS=`echo $3 | tr "+" " "`
		echo "$EXE $DIR/$2 $OPTS" >> $TMPNAME
		dos2unix $TMPNAME >> "$LOGS/cmds.txt" 2>&1;
		STATEMENT=`qsub $TMPNAME`;
		echo $STATEMENT >> "$LOGS/cmds.txt" 2>&1;
		HOLD=`echo $STATEMENT | cut -f3 -d" " |  cut -f1 -d"."`;
		echo $HOLD
	}
	
	####
	# Get Genes! <Must always be here>
	GENEHOLD=$(sendtoSGE $perl myGeneFinder.pl $DB/ucsc_refflat_hg19_2011-01-24.bed+$TMPDIR/VARFILE_\$SGE_TASK_ID 2G 1);
		
	###
	# Classify if they are known pathogenic variants (Boolean)
	P1_HOLD=$(sendtoSGE $perl pharmicoBool.pl $TMPDIR/VARFILE_\$SGE_TASK_ID+$DB/rsidSlim.tsv 2.5G 1);
	
	###
	# Query for HGMD Variant Information
	H1_HOLD=$(sendtoSGE $perl hgmdQuery_var.pl $TMPDIR/VARFILE_\$SGE_TASK_ID 512M 1);
	
	###
	# Query for Cosmic Variant Information
	C1_HOLD=$(sendtoSGE $perl cosmicFreq.pl $TMPDIR/VARFILE_\$SGE_TASK_ID 1G 1);
	
	
	##############################################
	# Gene Based Annotation -- Holds for Gene Annot Files
	#############################################
	#echo $GENEHOLD
	
	###
	# Query for Gene Test Information
	G2_HOLD=$(sendtoSGE $perl geneTestQuery.pl $TMPDIR/VARFILE_\${SGE_TASK_ID}_GENE 1G 1 $GENEHOLD);
		
	###
	# Query for Gene Information in HGMD
	H2_HOLD=$(sendtoSGE $perl hgmdQuery_gene.pl $TMPDIR/VARFILE_\${SGE_TASK_ID}_GENE 1G 1 $GENEHOLD);
	
	###
	# Query for Gene Information in Beauty Drug-Gene
	B2_HOLD=$(sendtoSGE java DrugTargetSingle_fat.jar -i+$TMPDIR/VARFILE_\${SGE_TASK_ID}_GENE+-o+$TMPDIR/VARFILE_\${SGE_TASK_ID}_BEAUTY+-db+$DB/druggable.sqlite 3G 1 $GENEHOLD);
	
	###
	# Query for Drugs, within a pathway of desired Gene
	DP2_HOLD=$(sendtoSGE java DrugTargetPathway.jar -i+$TMPDIR/VARFILE_\${SGE_TASK_ID}_GENE+-o+$TMPDIR/VARFILE_\${SGE_TASK_ID}_PATHWAY+-db+$DB/druggable2.sqlite 3G 1 $GENEHOLD);


	###
	# Query for CPIC Information, from PharmGKB
	C2_HOLD=$(sendtoSGE $perl cpicQuery.pl $TMPDIR/VARFILE_\${SGE_TASK_ID}_GENE 1G 1 $GENEHOLD);

	###
	# Query for Druggable Genome Information, from Sophic
	DG2_HOLD=$(sendtoSGE $perl druggable_flag.pl $TMPDIR/VARFILE_\${SGE_TASK_ID}_GENE+$DIR/DruggableGenome.txt 1G 1 $GENEHOLD);

	
	
	#######################################
	####  MERGE ANNOTATION COLLECTION  ####
	#######################################
	TABBEDOUT=$PREPATH"/"$RAWNAME"_Annotation.tsv";
	ALLHOLD="$P1_HOLD,$H1_HOLD,$C1_HOLD,$G2_HOLD,$H2_HOLD,$B2_HOLD,$C2_HOLD,$DG2_HOLD,$DP2_HOLD";
	
	###
	# Tool to merge intermediates back together, into a single file
	MERGE=$(sendtoSGE java FileJoiner.jar -d+$TMPDIR+-p+VARFILE_+-c+$NUM+-o+$TABBEDOUT 2G 0 $ALLHOLD);
		
	## Need to clean up ALL temp files!
	qsub -V -cwd -o $output_dir/logs -e $output_dir/logs -q $queue -m a -M $email -l h_stack=10M -l h_vmem=512M -hold_jid $MERGE $DIR/cleanUp.sh $TMPDIR
	
	
fi
