## March 20 2012

Problem: if there are less number of reads to realign for any Chromosome then GATK fails most of the time. 
Solution: workflow copies aligned bam as realigned bam with a warning message in the log






## Saurabh Baheti
## Feb 28th 2012

To run the workflow user needs to create 3 Configuration files i.e
1) Run information file
2) Tool information file
3) Sample information file
If it is a standard run then user should not change tool info and use the tool info in this space as it is. Make run info and sample info according to the samples you are running.

### how to run the workflow
After making these three config files user need to run this script to submit all the jobs

$ /projects/bsi/bictools/scripts/dev/WholeGenome/whole_genome.sh 
Usage: <Please specify path to run_info.txt file> 

$ qsub -V -cwd -q ngs-sec -l medp=TRUE -m a -M baheti.saurabh@mayo.edu /path/to/scripts/whole_genome.sh /path/to/run_info file

### after completion
After running the workflow, on completion user will reecieve an email about the completion of workflow, first thing user need to check is the errorlog and warning log file created on the output folder level for potenital errors.

#### after user validate and all teh data is correct and in right format then user can delete the intermediate files uisng this script

$ /path/to/scripts/cleanspace.sh /path/to/output_folder exome/whole_genom
$ /projects/bsi/bictools/scripts/dev/WholeGenome/cleanspace.sh 
Usage: </path/to/output dir > <tool (exome/whole_genome)>

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

EMIT_ALL_SITES=NO/YES YES if you want to output all the positions in the pileup 
if above paramerter is YES then specify these paramters
DEPTH_FILTER=4
TARGETTED=YES
PROB_FILTER=0.8



###########################################
