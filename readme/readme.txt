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
	drwxrwx--- 3 m088341 biostat  410 Sep  4 16:10 bam
drwxrws--- 2 m078940 biostat  148 Sep 16 20:46 input_fastq
-rwxrwx--- 1 tu03325 biostat  661 Sep 17 10:23 run_info.txt
-rwxrwx--- 1 tu03325 biostat  108 Sep 17 10:23 sample_info.txt
-rwxrwx--- 1 tu03325 biostat 7508 Sep 17 10:23 tool_info.txt
drwxrwx--- 2 m088341 biostat  312 Jan 31  2012 variant


##### For testing uisng bams
location of a dataset
$ /data2/bsi/RandD/WG_test/random/test_data
location of config files
$ /data2/bsi/RandD/WG_test/random/

#####OPTIONS for run info file

TOOL=GENOME_GPS (do not change this)
VERSION=1.0 (do not change this)
TYPE=exome/whole_genome
DISEASE=disease type (free text)
DATE=1/1/1(date of the analysis
READLENGTH=100 (should be a number)
PAIRED=1 (1 for paired end data and 0 for Single read data)
ANALYSIS=variant/external/mayo/annotation/alignment (For Sequnecing core run specify mayo to update secondary dashboard)
PI=baheti_saurabh (lastname_firstname_lanid)
MULTISAMPLE=NO (NO for single sample ; YES for paired analysis)
INPUT_DIR=/data2/bsi/RandD/WG_test/random/ (/path/to/input directory)
BASE_OUTPUT_DIR=/data2/bsi/RandD/WG_test/random/ (/path/to/output dirrectory) 
EMAIL=rstngsworkflow@mayo.edu (you should use only this email as we want to track the compute power and make html page using Asif's script )
USER_EMAIL=your email
SAMPLENAMES=s_normal:s_tumor (sample names ':' seperated)
GROUPNAMES=NA (if it is a paired analysis then provide the group names ':' seperated
LABINDEXES=-:- (index ':' seperated if user doesn't have index information then just put '-')
LANEINDEX=1:2 (lane information of a sample ':' seperated) [Keep in mind that #of samples=#of lanes=#od indexes]
CHRINDEX=1:2:3:4:5:6:7:8:9:10:11:12:13:14:15:16:17:18:19:20:21:22:X:Y:M (chr string user wants to use)
TOOL_INFO=/data2/bsi/RandD/WG_test/random/tool_info_all.txt (/path/to/tool_info file)
SAMPLE_INFO=/data2/bsi/RandD/WG_test/random/sample_info.bam.txt (/path/to/sample_info file)
OUTPUT_FOLDER=allModeule(output folder name i.e. flowcell id)
QUEUE=ngs-rand (queue information, ngs-sec -l medp ; ng-ext -l lowp)
CENTER=MAYO (do not change this)
PLATFORM=ILLUMINA (do not change this)
GENOMEBUILD=hg19(do not change this)
ALIGNER=NOVOALIGN/BWA
MARKDUP=YES/NO
FASTQC=YES/NO
FOLDER_FASTQC=NA or /path/to/fastqc results/
UPLOAD_TABLEBROWSER=NO/YES
REORDERSAM=NO/YES
VARIANT_TYPE=BOTH/SNV/INDEL
SNV_CALLER=GATK/SNVMIX
SOAMTIC_CALLER=SOMATICSNIPER/JOINTSNVMIX/MUTECT
SAMPLEINFORMATION=test run (free text )
DELIVERY_FOLDER=NA or /path/to/delivery folder to transfer realigned bams 
TERTIARY_FOLDER=NA or /path/to/tertiary folder for Biostats

### TOOL info file
### paramters with you can play with (user must provide these default values if there is no change)
MAX_READS_REALIGN=50000(default)
MAX_READS_MEM_REALIGN=100000(default) 
CNVNATOR_BINSIZE=1000(default)
PCT_READS_SEGSEQ=0.05(deafult)
MINFOLD=0.5(deafult)
MAXFOLD=1.5(deafult)
DISTGAP=1000(deafult)
BLAT_PORT=50000(deafult)
BLAT_SERVER=crick7(deafult)
HTTP_SERVER=bmidev2(deafult)
THREADS=4(deafult)
STRUCT_DIST_GENE=1000(deafult)
STRUCT_MIN_SUPPORT=10(default)
STRUCT_MIN_INDENTITY=0.9
T_DEPTH_FILTER=6
SOMATIC_THRESHOLD=0.1


EMIT_ALL_SITES=NO/YES YES if you want to output all the positions in the pileup 
if above paramerter is YES then specify these paramters
DEPTH_FILTER=4
TARGETTED=YES
PROB_FILTER=0.8


###########################################
