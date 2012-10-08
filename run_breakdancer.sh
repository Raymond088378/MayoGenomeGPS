#!/bin/bash

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
    echo -e "script to run break dancer for single or multiple samples\nUsage: samplenames input_bams </path/to/output directory> </path/to/run_info.txt> <group name>";
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
    script_path=$( cat $tool_info | grep -w '^WORKFLOW_PATH' | cut -d '=' -f2 )
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
	
########################################################	
######		

    export PERL5LIB=$perllib:$breakdancer;
    export PATH=$samtools:$PATH;
    cd $breakdancer

    if [ $multi != "YES" ]
    then
		mkdir -p $output_dir/$samples
		if (($SGE_TASK_ID <= $numchrs))
		then
			input_bam=$input/$samples/chr${chr}.cleaned.bam
			$samtools/samtools view -H $input_bam 1>$input_bam.break.header 2>$input_bam.fix.break.log
			if [ `cat $input_bam.fix.break.log | wc -l` -gt 0 ]
			then
				$script_path/email.sh $input_bam "bam is truncated or corrupt" $run_info
				$script_path/wait.sh $input_bam.fix.break.log 
			else
				rm $input_bam.fix.break.log
			fi	
			rm $input_bam.break.header
			out=$output_dir/$samples/
			## extarcting bam for $sample
			$samtools/samtools view -h $input_bam  | cut -f 1-11 | $samtools/samtools view -bS - > $out/$samples.$chr.bam
			$perl $breakdancer/bam2cfg.pl $out/$samples.$chr.bam > $out/$samples.$chr.cfg
			if [ -s $out/$samples.$chr.cfg ]
			then
				$breakdancer/breakdancer_max -c 5 -r 10 $out/$samples.$chr.cfg > $out/$samples.$chr.break
				cat $out/$samples.$chr.break |  sort -n -k 2,12n > $out/$samples.$chr.break.sorted
				#mv $out/$sample.$chr.break.sorted $out/$sample.$chr.break
				if [ -s $out/$samples.$chr.break.sorted ]
				then
					cat $out/$samples.$chr.break.sorted | awk '{print $1"\t"$2"\t"$2+1"\t"$4"\t"$5"\t"$5+1}' | $bedtools/pairToBed -a stdin -b $blacklist_sv -type neither | $script_path/report_original.pl $out/$samples.$chr.break.sorted > $out/$samples.$chr.break
				else
					touch $out/$samples.$chr.break
				fi	
				$script_path/Breakdancer2VCF.pl -i $out/$samples.$chr.break -f $ref_genome -o $out/$samples.$chr.break.vcf -s $samples -t $samtools
				$script_path/vcfsort.pl ${ref_genome}.fai $out/$samples.$chr.break.vcf > $out/$samples.$chr.break.vcf.sort
				mv $out/$samples.$chr.break.vcf.sort $out/$samples.$chr.break.vcf
				rm $out/$samples.$chr.cfg $out/$samples.$chr.bam

				if [ ! -s $out/$samples.$chr.break.vcf.fail ]
				then
					rm $out/$samples.$chr.break.vcf.fail
				fi  
			else
				$script_path/errorlog.sh $out/$samples.$chr.cfg run_cnvnator.sh ERROR "does not exist"
				touch $out/$samples.$chr.break
				touch $out/$samples.$chr.break.sorted
				$script_path/Breakdancer2VCF.pl -i $out/$samples.$chr.break -f $ref_genome -o $out/$samples.$chr.break.vcf -s $samples -t $samtools
				rm $out/$samples.$chr.break.vcf.fail
				$script_path/vcfsort.pl ${ref_genome}.fai $out/$samples.$chr.break.vcf > $out/$samples.$chr.break.vcf.sort
				mv $out/$samples.$chr.break.vcf.sort $out/$samples.$chr.break.vcf
			fi
		else	
			out=$output_dir/$samples/
			input_bam=$input/$samples.igv-sorted.bam
			$samtools/samtools view -H $input_bam 1>$input_bam.break.header 2>$input_bam.fix.break.log
			if [ `cat $input_bam.fix.break.log | wc -l` -gt 0 ]
			then
				$script_path/email.sh $input_bam "bam is truncated or corrupt" $run_info
				$script_path/wait.sh $input_bam.fix.break.log
			else
				rm $input_bam.fix.break.log
			fi	
			rm $input_bam.break.header
			$samtools/samtools view -h $input_bam | head -n 10000 | cut -f 1-11 | $samtools/samtools view -bS - > $out/$samples.tmp.bam
			if [ -s $out/$samples.tmp.bam ]
			then
				$perl $breakdancer/bam2cfg.pl $out/$samples.tmp.bam > $out/$samples.inter.cfg
				if [ -s $out/$samples.inter.cfg ]
				then
					rm $out/$samples.tmp.bam
					ln -s $input_bam $out/$samples.tmp.bam
					$breakdancer/breakdancer_max -t $out/$samples.inter.cfg > $out/$samples.inter.break
					cat $out/$samples.inter.break |  sort -n -k 2,12n > $out/$samples.inter.break.sorted
					# mv $out/$sample.inter.break.sorted $out/$sample.inter.break
					if [ -s $out/$samples.inter.break.sorted  ]
					then
						cat $out/$samples.inter.break.sorted | awk '{ print $1"\t"$2"\t"$2+1"\t"$4"\t"$5"\t"$5+1}' | $bedtools/pairToBed -a stdin -b $blacklist_sv -type neither | $script_path/report_original.pl $out/$samples.inter.break.sorted > $out/$samples.inter.break
					else
						touch $out/$samples.inter.break
					fi	
					#filter_sv_break $out/$samples.inter.break.sorted $out/$samples.inter.break
					$script_path/Breakdancer2VCF.pl -i $out/$samples.inter.break -f $ref_genome -o $out/$samples.inter.break.vcf -s $samples -t $samtools
					$script_path/vcfsort.pl ${ref_genome}.fai $out/$samples.inter.break.vcf > $out/$samples.inter.break.vcf.sort
					mv $out/$samples.inter.break.vcf.sort $out/$samples.inter.break.vcf
					if [ ! -s $out/$samples.inter.break.vcf.fail ]
					then
						rm $out/$samples.inter.break.vcf.fail
					fi  
					rm $out/$samples.inter.cfg $out/$samples.tmp.bam
				else
					$script_path/errorlog.sh $out/$samples.inter.cfg run_cnvnator.sh ERROR "does not exist"
					touch $out/$samples.inter.break
					touch $out/$samples.inter.break.sorted
					$script_path/Breakdancer2VCF.pl -i $out/$samples.inter.break -f $ref_genome -o $out/$samples.inter.break.vcf -s $samples -t $samtools
					rm $out/$samples.inter.break.vcf.fail
					$script_path/vcfsort.pl ${ref_genome}.fai $out/$samples.inter.break.vcf > $out/$samples.inter.break.vcf.sort
					mv $out/$samples.inter.break.vcf.sort $out/$samples.inter.break.vcf
				fi
			else
				$script_path/errorlog.sh $output_dir/$samples/$samples.tmp.bam run_cnvnator.sh ERROR "not created"
				exit 1;
			fi
		fi
	else
	    sample=$samples
	    mkdir -p $output_dir/$sample
	    input_bam=$input/$group.$sample.chr${chr}.bam
		$samtools/samtools view -H $input_bam 1>$input_bam.break.header 2>$input_bam.fix.break.log
		if [ `cat $input_bam.fix.break.log | wc -l` -gt 0 ]
		then
			$script_path/email.sh $input_bam "bam is truncated or corrupt" $run_info
			$script_path/wait.sh $input_bam.fix.break.log
		else
			rm $input_bam.fix.break.log
		fi	
		rm $input_bam.break.header
	    if (($SGE_TASK_ID <= $numchrs))
	    then	
            out=$output_dir/$sample/
            $samtools/samtools view -h $input_bam  | cut -f 1-11 | $samtools/samtools view -bS - >  $out/$sample.${chr}.bam
            $perl $breakdancer/bam2cfg.pl $out/$sample.$chr.bam > $out/$sample.$chr.cfg
		
            if [ -s $out/$sample.$chr.cfg ]
            then
                $breakdancer/breakdancer_max -c 5 -r 10 $out/$sample.$chr.cfg > $out/$sample.$chr.break
                cat $out/$sample.$chr.break |  sort -n -k 2,12n > $out/$sample.$chr.break.sorted
                if [ -s $out/$sample.$chr.break.sorted  ]
				then
					cat $out/$sample.$chr.break.sorted | awk '{ print $1"\t"$2"\t"$2+1"\t"$4"\t"$5"\t"$5+1}' | $bedtools/pairToBed -a stdin -b $blacklist_sv -type neither | $script_path/report_original.pl $out/$sample.$chr.break.sorted > $out/$sample.$chr.break
                else
					touch $out/$sample.$chr.break 
				fi
				$script_path/Breakdancer2VCF.pl -i $out/$sample.$chr.break -f $ref_genome -o $out/$sample.$chr.break.vcf -s $sample -t $samtools
                $script_path/vcfsort.pl ${ref_genome}.fai $out/$sample.$chr.break.vcf > $out/$sample.$chr.break.vcf.sort
                mv $out/$sample.$chr.break.vcf.sort $out/$sample.$chr.break.vcf
                rm $out/$sample.$chr.cfg $out/$sample.${chr}.bam
                if [ ! -s $out/$sample.$chr.break.vcf.fail ]
                then
                    rm $out/$sample.$chr.break.vcf.fail
                fi  
            else
                $script_path/errorlog.sh $out/$sample.$chr.cfg run_cnvnator.sh ERROR "not created"
                touch $out/$sample.$chr.break
				touch $out/$sample.$chr.break.sorted
				perl $script_path/Breakdancer2VCF.pl -i $out/$sample.$chr.break -f $ref_genome -o $out/$sample.$chr.break.vcf -s $sample -t $samtools
				rm $out/$sample.$chr.break.vcf.fail
                perl $script_path/vcfsort.pl ${ref_genome}.fai $out/$sample.$chr.break.vcf > $out/$sample.$chr.break.vcf.sort
                mv $out/$sample.$chr.break.vcf.sort $out/$sample.$chr.break.vcf
            fi
	    else	
            out=$output_dir/$sample/
            input_bam=$input/$group/${sample}.igv-sorted.bam
            $samtools/samtools view -H $input_bam 1>$input_bam.break.header 2>$input_bam.fix.break.log
			if [ `cat $input_bam.fix.break.log | wc -l` -gt 0 ]
			then
				$script_path/email.sh $input_bam "bam is truncated or corrupt" $run_info
				$script_path/wait.sh $input_bam.fix.break.log
			else
				rm $input_bam.fix.break.log
			fi
			rm $input_bam.break.header	
			$samtools/samtools view -h $input_bam | head -n 10000 | cut -f 1-11 | $samtools/samtools view -bS - > $out/$sample.tmp.bam
            if [ -s $out/$sample.tmp.bam ]
            then
                $perl $breakdancer/bam2cfg.pl $out/$sample.tmp.bam > $out/$sample.inter.cfg
                if [ -s $output_dir/$sample/$sample.inter.cfg ]
                then
                    rm $out/$sample.tmp.bam
                    ln -s $input_bam $out/$sample.tmp.bam   
                    $breakdancer/breakdancer_max -r 10 -t $out/$sample.inter.cfg > $out/$sample.inter.break
                    cat $out/$sample.inter.break |  sort -n -k 2,12n > $out/$sample.inter.break.sorted
                    if [ -s $out/$sample.inter.break.sorted ]
					then
						cat $out/$sample.inter.break.sorted | awk '{ print $1"\t"$2"\t"$2+1"\t"$4"\t"$5"\t"$5+1}' | $bedtools/pairToBed -a stdin -b $blacklist_sv -type neither | $script_path/report_original.pl $out/$sample.inter.break.sorted > $out/$sample.inter.break
                    else
						touch $out/$sample.inter.break
					fi
					perl $script_path/Breakdancer2VCF.pl -i $out/$sample.inter.break -f $ref_genome -o $out/$sample.inter.break.vcf -s $sample -t $samtools
                    perl $script_path/vcfsort.pl ${ref_genome}.fai $out/$sample.inter.break.vcf > $out/$sample.inter.break.vcf.sort
                    mv $out/$sample.inter.break.vcf.sort $out/$sample.inter.break.vcf
                    if [ ! -s $out/$sample.inter.break.vcf.fail ]
                    then
                        rm $out/$sample.inter.break.vcf.fail
                    fi  
                    rm $out/$sample.inter.cfg $out/$sample.tmp.bam
                else
					$script_path/errorlog.sh $out/$sample.inter.cfg run_cnvnator.sh ERROR "not created"
                    touch $out/$sample.inter.break
					touch $out/$sample.inter.break.sorted
					perl $script_path/Breakdancer2VCF.pl -i $out/$sample.inter.break -f $ref_genome -o $out/$sample.inter.break.vcf -s $sample -t $samtools
					rm  $out/$sample.inter.break.vcf.fail
                    perl $script_path/vcfsort.pl ${ref_genome}.fai $out/$sample.inter.break.vcf > $out/$sample.inter.break.vcf.sort
                    mv $out/$sample.inter.break.vcf.sort $out/$sample.inter.break.vcf
                fi
            else
                $script_path/errorlog.sh $out/$sample.tmp.bam run_cnvnator.sh ERROR "not created"
                exit 1;
            fi
	    fi
    fi
	echo `date`
fi
