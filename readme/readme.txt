### Saurabh Baheti
### Sept 17 2012

To run the workflow user needs to create 3 Configuration files i.e
1) Run information file
2) Tool information file
3) Sample information file
If it is a standard run then user should not change tool info and use the tool info in script folder as it is. Make run info and sample info according to the samples you are running.

### how to run the workflow
After making these three config files user need to run this script to submit all the jobs
Step 1. Generate teh unique identification number for a workflow instance
  	$ /projects/bsi/bictools/scripts/dnaseq/GENOME_GPS/trunk/unique_id.sh $run_info
  
Step 2. Run the workflow wrapper to submit jobs to the RCF cluster
	$ qsub -V -cwd -q ngs-sec -l medp=TRUE -m a -M baheti.saurabh@mayo.edu /path/to/scripts/whole_genome.sh /path/to/run_info_file
	$ /projects/bsi/bictools/scripts/dnaseq/GENOME_GPS/trunk/whole_genome.sh 
	Usage: <Please specify path to run_info.txt file> 

### after completion
After running the workflow, on completion user will reecieve an email about the completion of workflow, first thing user need to check is the errorlog and warning log file created on the output folder level for potenital errors.

#### after user validate and all the data is correct and in right format then user can delete the intermediate files using this script

	$ /path/to/scripts/transfer_clean.sh /path/to/output_folder /path/to/run_info_file
	$ /projects/bsi/bictools/scripts/dnaseq/GENOME_GPS/trunk/transfer_clean.sh
	Usage: wrapper to clean intermediate files and tansfer the data to tertiary, delivery folder
 	<secondary folder> < run_info >

### To test the workflow to get sample results user can use these test datasets
DIR : /data2/bsi/RandD/sampleData/treat 	
										/bam				:BAM files
										/input_fastq		:FASTQ files
										/variant			:VCF files

#####OPTIONS for run info file

TOOL=GENOME_GPS (do not change this)
VERSION=1.2 (do not change this)
TYPE=exome/whole_genome
DISEASE=disease type (free text)
READLENGTH=100 (should be a number)
PAIRED=1 (1 for paired end data and 0 for Single read data)
ANALYSIS=variant/external/mayo/annotation/alignment (For Sequnecing core run specify mayo to update secondary dashboard)
PI=baheti_saurabh_m078940 (lastname_firstname_lanid)
MULTISAMPLE=NO (NO for single sample ; YES for paired analysis)
INPUT_DIR=/data2/bsi/RandD/WG_test/random/ (/path/to/input directory)
BASE_OUTPUT_DIR=/data2/bsi/RandD/WG_test/random/ (/path/to/output dirrectory) 
SAMPLENAMES=s_normal:s_tumor (sample names ':' seperated)
GROUPNAMES=NA (if it is a paired analysis then provide the group names ':' seperated
LABINDEXES=-:- (index ':' seperated if user doesn't have index information then just put '-')
LANEINDEX=1:2 (lane information of a sample ':' seperated) [Keep in mind that #of samples=#of lanes=#od indexes]
CHRINDEX=1:2:3:4:5:6:7:8:9:10:11:12:13:14:15:16:17:18:19:20:21:22:X:Y:M (chr string user wants to use)
TOOL_INFO=/data2/bsi/RandD/WG_test/random/tool_info_all.txt (/path/to/tool_info file)
SAMPLE_INFO=/data2/bsi/RandD/WG_test/random/sample_info.bam.txt (/path/to/sample_info file)
OUTPUT_FOLDER=allModeule(output folder name i.e. flowcell id)
QUEUE=ngs-rand (queue information, ngs-sec -l medp ; ng-ext -l lowp)
GENOMEBUILD=hg19(do not change this)
ALIGNER=NOVOALIGN/BWA
FASTQC=YES/NO
FOLDER_FASTQC=NA or /path/to/fastqc results/
VARIANT_TYPE=BOTH/SNV/INDEL
SNV_CALLER=GATK/SNVMIX/BEAUTY_EXOME
SOAMTIC_CALLER=SOMATICSNIPER/JOINTSNVMIX/MUTECT/BEAUTY_EXOME
SAMPLEINFORMATION=test run (free text )
DELIVERY_FOLDER=NA or /path/to/delivery folder to transfer realigned bams for visualization
TERTIARY_FOLDER=NA or /path/to/tertiary folder for Biostats

### TOOL info file
### paramters with you can play with (user must provide these default values if there is no change)
THREADS=4
REORDERSAM=NO
EMIT_ALL_SITES=NO
VARIANT_FILTER=YES
SOMATIC_VARIANT_FILTER=YES
DEPTH_FILTER=0
TARGETTED=YES
MARKDUP=YES
REMOVE_DUP=FALSE
REMOVE_ALIGNED_BAM=YES
T_DEPTH_FILTER=6
UPLOAD_TABLEBROWSER=YES
PLATFORM=illumina
CENTER=mayo
QUEUE=ngs-rand

#Parameters
SNVMIX2_params=
SNVMIX2_Filter=-p 0.8
UnifiedGenotyper_params=-maxAlleles 5
SOMATIC_INDEL_params=--window_size 1000
SOMATIC_SNIPER_params=-q 20 -Q 20
MUTECT_params=
BREAKDANCER_params=-c 5 -r 10
CREST_params=
JSM_Filter=-prob 0.1
JOINTSNVMIX_params=
NOVO_params=-g 60 -x 2 -i PE 425,80 -r Random --hdrhd off -v 120
BWA_params=-l 32 -t 4
VQSR_params_SNV=--maxGaussians 4 --percentBadVariants 0.05
VQSR_params_INDEL=--maxGaussians 4 --percentBadVariants 0.12
PICARD_ReadGroup_params=PL=illumina CN=mayo LB=hg19 CREATE_INDEX=true


MAX_READS_REALIGN=50000
MAX_READS_MEM_REALIGN=150000
CNVNATOR_BINSIZE=1000
PCT_READS_SEGSEQ=0.05
MINFOLD=0.5
MAXFOLD=1.5
DISTGAP=1000
BLAT_PORT=50000
BLAT_SERVER=localhost
STRUCT_DIST_GENE=1000
STRUCT_MIN_SUPPORT=10
STRUCT_MIN_IDENTITY=0.9
STRUCT_PCT_BLACKLIST=1
SNP_DISTANCE_INDEL=10
MAX_FILE_HANDLES=100
MAX_READS_MEM_SORT=2000000
WINDOW_BLAT=50
TB_PORT=8886
TB_HOST=charlotte
JOB_LIMIT=3000


###########################################
