#!/bin/sh

########################################################
###### 	SV CALLER FOR WHOLE GENOME ANALYSIS PIPELINE

######		Program:			run_cnvnator.sh
######		Date:				09/26/2011
######		Summary:			Calls CNVnator
######		Input 
######		$1	=	samplename
######		$2	=	input_bam
######		$3	=	/path/to/output directory
######		$4	=	/path/to/run_info.txt
######		Output files:	BAM files. 
#####		added -r 10 parmeter min supported reads
########################################################

if [ $# -le 3 ]
then
    echo -e "Usage: script to run break dancer for single or multiple samples \n samplenames input_bams </path/to/output directory> </path/to/run_info.txt> <group name>";
else
    set -x
    echo `date`
    samples=$1
    input=$2
    output_dir=$3
    run_info=$4
    if [ $5 ]
    then
        group=$5
    fi
#SGE_TASK_ID=1

########################################################	
######		Reading run_info.txt and assigning to variables

    tool_info=$( cat $run_info | grep -w '^TOOL_INFO' | cut -d '=' -f2)
    sample_info=$( cat $run_info | grep -w '^SAMPLE_INFO' | cut -d '=' -f2)
    script_path=$( cat $tool_info | grep -w '^WHOLEGENOME_PATH' | cut -d '=' -f2 )
    breakdancer=$( cat $tool_info | grep -w '^BREAKDANCER' | cut -d '=' -f2 )
    perl=$( cat $tool_info | grep -w '^PERL_BREAKDANCER' | cut -d '=' -f2 )
    perllib=$( cat $tool_info | grep -w '^PERLLIB_BREAKDANCER' | cut -d '=' -f2 )
    samtools=$( cat $tool_info | grep -w '^SAMTOOLS' | cut -d '=' -f2 )
    numchrs=$(cat $run_info | grep -w '^CHRINDEX' | cut -d '=' -f2 | tr ":" "\n" | wc -l)
    chr=$(cat $run_info | grep -w '^CHRINDEX' | cut -d '=' -f2 | tr ":" "\n" | head -n $SGE_TASK_ID | tail -n 1)
    ref_genome=$( cat $tool_info | grep -w '^REF_GENOME' | cut -d '=' -f2 )
    multi=$( cat $run_info | grep -w '^MULTISAMPLE' | cut -d '=' -f2| tr "[a-z]" "[A-Z]")

    blacklist_sv=$( cat $tool_info | grep -w '^BLACKLIST_SV' | cut -d '=' -f2 )
    bedtools=$( cat $tool_info | grep -w '^BEDTOOLS' | cut -d '=' -f2 )
	PATH=$bedtools:$PATH
#filter_sv_break ()
#{
    # Removes SV from file ($1) that intersect $blacklist_sv and puts the result in
    # outfile ($2)
#    file=$1
#    outfile=$2

#    cat $file | awk '{ print $1"\t"$2"\t"$2+1"\t"$4"\t"$5"\t"$5+1}' |\ 
#    $bedtools/pairToBed -a stdin -b $blacklist_sv -type neither |\ 
#    $script_path/report_original.pl $file > $outfile
#}


########################################################	
######		

    echo `date`
    PERL5LIB=$perllib;
    PATH=$samtools:$PATH;

#i=1
#   for sample in `echo $samples | tr ":" "\n"`
#   do
#     sampleArray[$i]=$sample
#     let i=i+1
#   done

    if [ $multi != "YES" ]
    then
	mkdir -p $output_dir/$samples
	if (($SGE_TASK_ID <= $numchrs))
	then
	    input_bam=$input/$samples/chr${chr}.cleaned.bam
	    out=$output_dir/$samples/
	    ## extarcting bam for $sample
	    $samtools/samtools view -h $input_bam  | cut -f 1-11 | $samtools/samtools view -bS - > $out/$samples.$chr.bam
	    $perl $breakdancer/perl/bam2cfg.pl $out/$samples.$chr.bam > $out/$samples.$chr.cfg
			
	    if [ -s $out/$samples.$chr.cfg ]
	    then
		$breakdancer/cpp/breakdancer_max -c 5 -r 10 $out/$samples.$chr.cfg > $out/$samples.$chr.break
		cat $out/$samples.$chr.break |  sort -n -k 2,12n > $out/$samples.$chr.break.sorted
		#mv $out/$sample.$chr.break.sorted $out/$sample.$chr.break

		cat $out/$samples.$chr.break.sorted | awk '{print $1"\t"$2"\t"$2+1"\t"$4"\t"$5"\t"$5+1}' | $bedtools/pairToBed -a stdin -b $blacklist_sv -type neither | $script_path/report_original.pl $out/$samples.$chr.break.sorted > $out/$samples.$chr.break

		perl $script_path/Breakdancer2VCF.pl -i $out/$samples.$chr.break -f $ref_genome -o $out/$samples.$chr.break.vcf -s $samples -t $samtools
		perl $script_path/vcfsort.pl ${ref_genome}.fai $out/$samples.$chr.break.vcf > $out/$samples.$chr.break.vcf.sort
		mv $out/$samples.$chr.break.vcf.sort $out/$samples.$chr.break.vcf
		rm $out/$samples.$chr.cfg $out/$samples.$chr.bam

		if [ ! -s $out/$samples.$chr.break.vcf.fail ]
		then
		    rm $out/$samples.$chr.break.vcf.fail
		fi  
	    else
			echo "Error Breakdancer: File $out/$sample.$chr.cfg not created"
			exit 1
	    fi
    else	
		out=$output_dir/$samples/
		input_bam=$input/$samples.igv-sorted.bam
		$samtools/samtools view -h $input_bam | head -n 10000 | cut -f 1-11 | $samtools/samtools view -bS - > $out/$samples.tmp.bam
		if [ -s $out/$samples.tmp.bam ]
		then
			$perl $breakdancer/perl/bam2cfg.pl $out/$samples.tmp.bam > $out/$samples.inter.cfg
			if [ -s $out/$samples.inter.cfg ]
			then
				rm $out/$samples.tmp.bam
				ln -s $input_bam $out/$samples.tmp.bam
				$breakdancer/cpp/breakdancer_max -t $out/$samples.inter.cfg > $out/$samples.inter.break
				cat $out/$samples.inter.break |  sort -n -k 2,12n > $out/$samples.inter.break.sorted
				# mv $out/$sample.inter.break.sorted $out/$sample.inter.break
				cat $out/$samples.inter.break.sorted | awk '{ print $1"\t"$2"\t"$2+1"\t"$4"\t"$5"\t"$5+1}' | $bedtools/pairToBed -a stdin -b $blacklist_sv -type neither | $script_path/report_original.pl $out/$samples.inter.break.sorted > $out/$samples.inter.break
				#filter_sv_break $out/$samples.inter.break.sorted $out/$samples.inter.break
		
				perl $script_path/Breakdancer2VCF.pl -i $out/$samples.inter.break -f $ref_genome -o $out/$samples.inter.break.vcf -s $samples -t $samtool

				perl $script_path/vcfsort.pl ${ref_genome}.fai $out/$samples.$chr.break.vcf > $out/$samples.$chr.break.vcf.sort
				mv $out/$samples.$chr.break.vcf.sort $out/$samples.$chr.break.vcf
				if [ ! -s $out/$samples.inter.break.vcf.fail ]
				then
					rm $out/$samples.inter.break.vcf.fail
				fi  
				rm $out/$samples.inter.cfg $out/$samples.tmp.bam
			else
				echo "Error Breakdancer: File $output_dir/$samples/$samples.inter.cfg not created"
				exit 1
			fi
		else
			echo "Error Breakdancer: File $output_dir/$samples/$samples.tmp.bam not created" 
			exit 1
		fi
	fi
	else
#for i in `cat $sample_info| grep -w "^$samples" | cut -d '=' -f2`
#	  do
	    sample=$samples
	    mkdir -p $output_dir/$sample
	    input_bam=$input/$group.$sample.chr${chr}.bam
	    if (($SGE_TASK_ID <= $numchrs))
	    then	
            out=$output_dir/$sample/
            $samtools/samtools view -h $input_bam  | cut -f 1-11 | $samtools/samtools view -bS - >  $out/$sample.${chr}.bam
            $perl $breakdancer/perl/bam2cfg.pl $out/$sample.$chr.bam > $out/$sample.$chr.cfg
		
            if [ -s $out/$sample.$chr.cfg ]
            then
                $breakdancer/cpp/breakdancer_max -c 5 -r 10 $out/$sample.$chr.cfg > $out/$sample.$chr.break
                cat $out/$sample.$chr.break |  sort -n -k 2,12n > $out/$sample.$chr.break.sorted
                #mv $out/$sample.$chr.break.sorted $out/$sample.$chr.break
                # filter_sv_break $out/$sample.$chr.break.sorted $out/$sample.$chr.break
                cat $out/$sample.$chr.break.sorted | awk '{ print $1"\t"$2"\t"$2+1"\t"$4"\t"$5"\t"$5+1}' | $bedtools/pairToBed -a stdin -b $blacklist_sv -type neither | $script_path/report_original.pl $out/$sample.$chr.break.sorted > $out/$sample.$chr.break
                perl $script_path/Breakdancer2VCF.pl -i $out/$sample.$chr.break -f $ref_genome -o $out/$sample.$chr.break.vcf -s $sample -t $samtools
                perl $script_path/vcfsort.pl ${ref_genome}.fai $out/$sample.$chr.break.vcf > $out/$sample.$chr.break.vcf.sort
                mv $out/$sample.$chr.break.vcf.sort $out/$sample.$chr.break.vcf
                rm $out/$sample.$chr.cfg $out/$sample.${chr}.bam
                if [ ! -s $out/$sample.$chr.break.vcf.fail ]
                then
                    rm $out/$sample.$chr.break.vcf.fail
                fi  
            else
                echo "Error Breakdancer: File $out/$sample.$chr.cfg not created"
                exit 1
            fi
	    else	
            out=$output_dir/$sample/
            input_bam=$input/${sample}.igv-sorted.bam
            $samtools/samtools view -h $input_bam | head -n 10000 | cut -f 1-11 | $samtools/samtools view -bS - > $out/$sample.tmp.bam

            if [ -s $out/$sample.tmp.bam ]
            then
                $perl $breakdancer/perl/bam2cfg.pl $out/$sample.tmp.bam > $out/$sample.inter.cfg
                if [ -s $output_dir/$sample/$sample.inter.cfg ]
                then
                    rm $out/$sample.tmp.bam
                    ln -s $input_bam $out/$sample.tmp.bam   
                    $breakdancer/cpp/breakdancer_max -r 10 -t $out/$sample.inter.cfg > $out/$sample.inter.break
                    cat $out/$sample.inter.break |  sort -n -k 2,12n > $out/$sample.inter.break.sorted
                    # mv $out/$sample.inter.break.sorted $out/$sample.inter.break
                    #filter_sv_break $out/$sample.inter.break.sorted $out/$sample.inter.break
                    cat $out/$sample.inter.break.sorted | awk '{ print $1"\t"$2"\t"$2+1"\t"$4"\t"$5"\t"$5+1}' | $bedtools/pairToBed -a stdin -b $blacklist_sv -type neither | $script_path/report_original.pl $out/$sample.inter.break.sorted > $out/$sample.inter.break
                    perl $script_path/Breakdancer2VCF.pl -i $out/$sample.inter.break -f $ref_genome -o $out/$sample.inter.break.vcf -s $sample -t $samtools
                    perl $script_path/vcfsort.pl ${ref_genome}.fai $out/$sample.$chr.break.vcf > $out/$sample.$chr.break.vcf.sort
                    mv $out/$sample.$chr.break.vcf.sort $out/$sample.$chr.break.vcf
                    if [ ! -s $out/$sample.inter.break.vcf.fail ]
                    then
                        rm $out/$sample.inter.break.vcf.fail
                    fi  
                    rm $out/$sample.inter.cfg $out/$sample.tmp.bam
                else
                    echo "Error Breakdancer: File $output_dir/$sample/$sample.inter.cfg not created"
                    exit 1
                fi
            else
                echo "Error Breakdancer: File $output_dir/$sample/$sample.tmp.bam not created" 
                exit 1
            fi
	    fi
# done
    fi
fi
