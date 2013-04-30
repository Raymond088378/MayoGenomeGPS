#!/bin/bash
##	INFO
#	to cat all the realigned BAM for IGV visualization

######################################
#		$1		=	input folder (realignment sample folder)
#		$3		=	sample
#		$4		=	output folder
#		$5		=	run info file
#########################################

if [ $# -le 4 ];
then
    echo -e "SCRIPT to create IGV BAM \
		\nUsage: ./igv_bam.sh </path/to/realign dir> </path/to/output folder>  \
			<sample> </path/to/alignment folder></path/to/run info> <group name (optional)";
else	
    set -x
    echo `date`
    input=$1
    output=$2
    sample=$3
    alignment=$4
    run_info=$5
    
    if [ $6 ]
    then
    	group=$6
	fi
	
    tool_info=$( cat $run_info | grep -w '^TOOL_INFO' | cut -d '=' -f2)
    sample_info=$( cat $run_info | grep -w '^SAMPLE_INFO' | cut -d '=' -f2)
    chrs=$( cat $run_info | grep -w '^CHRINDEX' | cut -d '=' -f2)
    chrIndexes=$( echo $chrs | tr ":" "\n" )
    samtools=$( cat $tool_info | grep -w '^SAMTOOLS' | cut -d '=' -f2 )
    script_path=$( cat $tool_info | grep -w '^WORKFLOW_PATH' | cut -d '=' -f2 )	
	delivery_folder=$( cat $run_info | grep -w '^DELIVERY_FOLDER' | cut -d '=' -f2)
    multi=$( cat $run_info | grep -w '^MULTISAMPLE' | cut -d '=' -f2| tr "[a-z]" "[A-Z]")
	analysis=$( cat $run_info | grep -w '^ANALYSIS' | cut -d '=' -f2| tr "[A-Z]" "[a-z]" )
	remove_bam=$( cat $tool_info | grep -w '^REMOVE_ALIGNED_BAM' | cut -d '=' -f2 | tr "[a-z]" "[A-Z]")
    stop_after_realignment=$( cat $tool_info | grep -w '^STOP_AFTER_REALIGNMENT' | cut -d '=' -f2 | tr "[a-z]" "[A-Z]" )
	i=1
    for chr in $chrIndexes
    do
        chrArray[$i]=$chr
        let i=i+1
    done
	
	if [ $analysis == "variant" ]
	then
		previous="split_bam_chr.sh"
	else
		previous="realign_recal.sh"
	fi
	
    if [ $6 ]
    then
        cd $input/$group
		mkdir -p $output/$group
        samples=$( cat $sample_info | grep -w "^$group" | cut -d '=' -f2)
        for j in $(seq 1 ${#chrArray[@]})
		do
			chr=${chrArray[$j]}
			i=chr$chr.cleaned.bam
			$samtools/samtools view -H $i 1>$i.igv.header 2> $i.igv.fix.log
			if [[ `cat $i.igv.fix.log | wc -l` -gt 0 || `cat $i.igv.header | wc -l` -le 0  ]]
			then
				$script_path/email.sh $i "bam is truncated or corrupt" $previous $run_info
				$script_path/wait.sh $i.igv.fix.log
			else
				rm $i.igv.fix.log
			fi	
			rm $i.igv.header
        done

    	sam=`echo $samples | tr " " "\n" | grep -w -v $sample | tr "\n" " "`
        gr=""
    	for s in $sam
    	do
        	a="ID:$s|";
        	gr="$gr$a"
    	done
    	gr=`echo $gr |  sed "s/|$//"`
    	cat $output/$group.$sample.chr$chr.bam
    	cat $output/$group/$sample.header.sam |grep -w -E -v "$gr" > $output/$sample/$sample.$i.header.sam
            
            if [ ${#chrArray[@]} -gt 1 ]
            then
                input_bam=""
                for j in $(seq 1 ${#chrArray[@]})
                do
                    chr=${chrArray[$j]}
                    if [ -f $output/$sample.$i.chr$chr.bam ]
					then
						input_bam="$input_bam $output/$sample.$i.chr$chr.bam"
					fi
				done
				$samtools/samtools merge -h $output/$sample/$sample.$i.header.sam $output/$sample/$i.igv-sorted.bam $input_bam $output/$i.sorted.bam.extra.bam 
			else
                chr=${chrArray[1]}
				$samtools/samtools view -H $output/$sample.$i.chr$chr.bam > $output/$sample/$sample.1.header.sam
				$samtools/samtools merge -h $output/$sample/$sample.1.header.sam $output/$sample/$i.igv-sorted.bam $output/$sample.$i.chr$chr.bam $output/$i.sorted.bam.extra.bam 
				rm $output/$sample/$sample.1.header.sam
            fi
            if [ -s $output/$sample/$i.igv-sorted.bam ]
            then
                $samtools/samtools index $output/$sample/$i.igv-sorted.bam
                if [ -s $output/$sample/$sample.$i.header.sam ]
				then
					rm $output/$sample/$sample.$i.header.sam
				fi
				if [ -s $output/$i.sorted.bam.extra.bam ]
				then
					rm $output/$i.sorted.bam.extra.bam $output/$i.sorted.bam.extra.bam.bai
				fi
				if [ $analysis != "variant" ]
				then
					if [ $remove_bam == "YES" ]
					then
						rm $alignment/$i/$i.sorted.bam $alignment/$i/$i.sorted.bam.bai
					fi
				fi
			else
                $script_path/errorlog.sh $output/$i.igv-sorted.bam igv_bam.sh ERROR "failed to create"
                exit 1;
            fi
        done
		rm $output/$sample/$sample.header.sam
		
    else
        cd $input/$sample/
        for j in $(seq 1 ${#chrArray[@]})
		do
			chr=${chrArray[$j]}
			i=chr$chr.cleaned.bam
            $samtools/samtools view -H $i 1>$i.igv.header 2> $i.fix.igv.log
			if [[ `cat $i.fix.igv.log | wc -l` -gt 0 || `cat $i.igv.header | wc -l` -le 0 ]]
			then
				$script_path/email.sh $i "bam is truncated or corrupt" $previous $run_info
				$script_path/wait.sh $i.fix.igv.log
			else
				rm $i.fix.igv.log
			fi
			rm $i.igv.header	
        done
        # only merge if there is more than 1 chr
        ### hard coding to find extension *.cleaned.bam
        if [ ${#chrArray[@]} -gt 1 ]
        then
            input_bam=""
            index=""
            for i in $(seq 1 ${#chrArray[@]})
            do
                if [ -f $input/$sample/chr${chrArray[$i]}.cleaned.bam ]
				then
					input_bam="$input_bam $input/$sample//chr${chrArray[$i]}.cleaned.bam"
					index="$index $input/$sample/chr${chrArray[$i]}.cleaned.bam.bai"
				fi
			done
			$script_path/sortbam.sh "$input_bam $output/$sample.sorted.bam.extra.bam" $output/$sample.igv-sorted.bam $output coordinate true $run_info no yes
        else
			$script_path/sortbam.sh "$input/$sample/chr${chrArray[1]}.cleaned.bam $output/$sample.sorted.bam.extra.bam" \
				$output/$sample.igv-sorted.bam $output coordinate true $run_info no yes
        fi
        
        if [ -s $output/$sample.igv-sorted.bam ]
        then
            $script_path/indexbam.sh $output/$sample.igv-sorted.bam $tool_info
            
            if [ -s $output/$sample.sorted.bam.extra.bam ]
			then
				rm $output/$sample.sorted.bam.extra.bam $output/$sample.sorted.bam.extra.bam.bai
			fi	
			if [ $analysis != "variant" ]
			then
				if [ $remove_bam == "YES" ]
				then
					rm $alignment/$sample/$sample.sorted.bam $alignment/$sample/$sample.sorted.bam.bai
				fi
			fi
		else
            $script_path/errorlog.sh $output/$sample.igv-sorted.bam igv_bam.sh ERROR "failed to create"
			exit 1;
        fi
    fi     
    out=$delivery_folder/IGV_BAM
    if [ $delivery_folder != "NA" ]
    then
        if [ -d $delivery_folder ]
        then
            if [ ! -d $out ]
            then
                mkdir -p $out
				chmod -Rf 777 $out
            fi
            if [ $multi == "YES" ]
            then
                mkdir -p $out/$sample
				chmod -Rf 777 $out/$sample
				pair=$( cat $sample_info | grep -w "$sample" | cut -d '=' -f2)
                for i in $pair
                do
                    mv $output/$sample/$i.igv-sorted.bam $out/$sample/
                    ln -s $out/$sample/$i.igv-sorted.bam $output/$sample/$i.igv-sorted.bam
                    mv $output/$sample/$i.igv-sorted.bam.bai $out/$sample/
                    ln -s $out/$sample/$i.igv-sorted.bam.bai $output/$sample/$i.igv-sorted.bam.bai
                done
            else
                mv $output/$sample.igv-sorted.bam $out/
                ln -s $out/$sample.igv-sorted.bam $output/$sample.igv-sorted.bam
                mv $output/$sample.igv-sorted.bam.bai $out/
                ln -s $out/$sample.igv-sorted.bam.bai $output/$sample.igv-sorted.bam.bai
            fi
        else
            $script_path/errorlog.sh $delivery_folder igv_bam.sh ERROR "folder not exist"
            exit 1
        fi
    fi  
	if [ $stop_after_realignment == "YES" ]
    then
		if [ $multi == "YES" ]
		then
			sam=$( cat $sample_info | grep -w "^$sample" | cut -d '=' -f2 | tr "\t" " ")
			for samp in $sam
			do
				$script_path/dashboard.sh $samp $run_info Complete complete
			done
		else
			$script_path/dashboard.sh $sample $run_info Complete complete
		fi
    fi
	echo `date`
fi
