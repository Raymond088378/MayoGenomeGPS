#!/bin/bash

#	INFO
#	the script generate a merged report per sample includes sift and sseq annotation (getting rid of redundant columns )

########################################
#		$1		=	Temp Reports folder
#		$2		=	sample name
#		$3		=	chromosome
#		$4		=	sseq directory
#		$5		=	sift directory
#		$6		=	variant file
#		$7		=	run info
############################################
		
if [ $# != 8 ];
then
	echo "Usage:<Tempreports><sample name><chromosome><sseq dir><sift dir><variant file><run info>";
else	
	set -x
	echo `date`
	TempReports=$1
	sample=$2 
	which_chr=$3
	sift=$4 
	snpeff=$5
	poly=$6
	var=$7	
	run_info=$8
	tool_info=$( cat $run_info | grep -w '^TOOL_INFO' | cut -d '=' -f2)
	script_path=$( cat $tool_info | grep -w '^WORKFLOW_PATH' | cut -d '=' -f2 )
	codon_ref=$( cat $tool_info | grep -w '^CODON_REF' | cut -d '=' -f2)
	GeneIdMap=$( cat $tool_info | grep -w '^GeneIdMap' | cut -d '=' -f2)
	UCSC=$( cat $tool_info | grep -w '^UCSC_TRACKS' | cut -d '=' -f2 )
	bed=$( cat $tool_info | grep -w '^BEDTOOLS' | cut -d '=' -f2 )
	blacklisted=$( cat $tool_info | grep -w '^BLACKLISTED' | cut -d '=' -f2 )
	mapability=$( cat $tool_info | grep -w '^MAPABILITY' | cut -d '=' -f2 )
	repeat=$( cat $tool_info | grep -w '^REPEATREGIONS' | cut -d '=' -f2 )
	miRbase=$( cat $tool_info | grep -w '^miRbase' | cut -d '=' -f2 )
	ssr=$( cat $tool_info | grep -w '^SNP_SR' | cut -d '=' -f2 )
	scs=$( cat $tool_info | grep -w '^SNP_CS' | cut -d '=' -f2 )
	sao=$( cat $tool_info | grep -w '^SNP_SAO' | cut -d '=' -f2 )
	build=$( cat $tool_info | grep -w '^SNP_BUILD' | cut -d '=' -f2 )
	typeset -i codon
	typeset -i SNP_Type
	
	### add sift columns
	file=$TempReports/$var.rsIDs.allele_frequencies
	num=`cat $file | awk '{print $1"_"$2}'| wc -l`
	chr=`awk -F '\t' '{ for(i=1;i<=NF;i++){ if ($i == "Chr") {print i} } }' $file`
	pos=`awk -F '\t' '{ for(i=1;i<=NF;i++){ if ($i == "Position") {print i} } }' $file`
	ref=`awk -F '\t' '{ for(i=1;i<=NF;i++){ if ($i == "Ref") {print i} } }' $file`
	alt=`awk -F '\t' '{ for(i=1;i<=NF;i++){ if ($i == "Alt") {print i} } }' $file`
	if [ -s $sift/${sample}_chr${which_chr}_predictions.tsv ]
	then
		perl $script_path/parse_siftPredictions.pl -i $file -s $sift/${sample}_chr${which_chr}_predictions.tsv -c $chr -p $pos -r $ref -a $alt -o $file.sift
	else
		$script_path/email.sh $sift/${sample}_chr${which_chr}_predictions.tsv  "not found" $JOB_NAME $JOB_ID $run_info
		touch $sift/${sample}_chr${which_chr}_predictions.tsv.fix.log
		while [ -f $sift/${sample}_chr${which_chr}_predictions.tsv.fix.log ]
		do
			echo "waiting for the fix"
			sleep 2m
		done
		perl $script_path/parse_siftPredictions.pl -i $file -s $sift/${sample}_chr${which_chr}_predictions.tsv -c $chr -p $pos -r $ref -a $alt -o $file.sift
	fi
	
	### sort the erport using the start and position
	perl $script_path/sort.variantReport.pl -i $file.sift -o $file.sift.sort -f Position
	mv $file.sift.sort $file.sift
	num_a=`cat $file.sift | wc -l`
	if [ $num == $num_a ]
	then
	    rm $file
	else
	    echo "ERROR : sift addition failed for $file"
	fi    
	##  add codon preference
	codon=`awk -F '\t' '{ for(i=1;i<=NF;i++){ if ($i ~ "Codons") {print i} } }' $file.sift`
	SNP_Type=`awk -F '\t' '{ for(i=1;i<=NF;i++){ if ($i ~ "SNP Type") {print i} } }' $file.sift`
	cat $file.sift | cut -f ${codon},${SNP_Type} > $file.sift.forCodons.2columns
	perl $script_path/codon_pref.pl $codon_ref $file.sift.forCodons.2columns $file.sift.forCodons.2columns.added
	paste $file.sift $file.sift.forCodons.2columns.added > $file.sift.codons
	num=`cat $file.sift.codons | wc -l`
	if [ $num_a == $num ]
	then
	    rm $file.sift
	    rm $file.sift.forCodons.2columns.added
	    rm $file.sift.forCodons.2columns
	else
	    echo "ERROR : adding codon prefernce failed for $file"
	fi    
	### blacklisted, Alignablility or Uniqueness columns
	cat $file.sift.codons | awk 'NR>2' | awk '{print $1"\t"($2-1)"\t"$2}' > $file.sift.codons.ChrPos.bed
	if [ -s $file.sift.codons.ChrPos.bed ]
	then
		$bed/intersectBed -a $file.sift.codons.ChrPos.bed -b $blacklisted -c | awk '{print $NF}' > $file.sift.codons.ChrPos.bed.black.i
	else
		touch $file.sift.codons.ChrPos.bed.black.i
	fi	
	echo -e "\nBlacklistedRegion" > $file.sift.codons.ChrPos.bed.black.i.tmp
	cat  $file.sift.codons.ChrPos.bed.black.i.tmp $file.sift.codons.ChrPos.bed.black.i > $file.sift.codons.ChrPos.bed.black
	rm $file.sift.codons.ChrPos.bed.black.i.tmp $file.sift.codons.ChrPos.bed.black.i
	if [ -s $file.sift.codons.ChrPos.bed ]
	then
		$bed/intersectBed -a $file.sift.codons.ChrPos.bed -b $mapability -c | awk '{print $NF}' > $file.sift.codons.ChrPos.bed.map.i
	else
		touch $file.sift.codons.ChrPos.bed.map.i
	fi
	echo -e "\nAlignability/Uniquness" > $file.sift.codons.ChrPos.bed.map.i.tmp
	cat  $file.sift.codons.ChrPos.bed.map.i.tmp $file.sift.codons.ChrPos.bed.map.i > $file.sift.codons.ChrPos.bed.map
	rm $file.sift.codons.ChrPos.bed.map.i.tmp $file.sift.codons.ChrPos.bed.map.i
	paste $file.sift.codons $file.sift.codons.ChrPos.bed.black $file.sift.codons.ChrPos.bed.map > $file.sift.codons.map
	rm $file.sift.codons.ChrPos.bed $file.sift.codons.ChrPos.bed.black $file.sift.codons.ChrPos.bed.map
	num_a=`cat $file.sift.codons.map | wc -l`
	if [ $num == $num_a ]
	then
	    rm $file.sift.codons
	else
	    echo "ERROR : adding blacklisted, Alignablility or Uniqueness columns failed for $file"
		exit 1;
	fi    
	## intersect with repeat region bed file
	cat $file.sift.codons.map | awk 'NR>2' | awk '{print $1"\t"($2-1)"\t"$2}' > $file.sift.codons.map.ChrPos.bed
	if [ -s $file.sift.codons.map.ChrPos.bed ]
	then
		$bed/intersectBed -a $file.sift.codons.map.ChrPos.bed -b $repeat -c | awk '{print $NF}' > $file.sift.codons.map.ChrPos.bed.i
	else
		touch $file.sift.codons.map.ChrPos.bed.i
	fi
	echo -e "\nRepeat_Region" >> $file.sift.codons.map.ChrPos.bed.i.tmp
	cat $file.sift.codons.map.ChrPos.bed.i >> $file.sift.codons.map.ChrPos.bed.i.tmp
	paste $file.sift.codons.map $file.sift.codons.map.ChrPos.bed.i.tmp > $file.sift.codons.map.repeat
	rm $file.sift.codons.map.ChrPos.bed.i.tmp $file.sift.codons.map.ChrPos.bed.i $file.sift.codons.map.ChrPos.bed
	num=`cat $file.sift.codons.map.repeat | wc -l`
	if [ $num_a == $num ]
	then
	    rm $file.sift.codons.map
	else
	    echo "ERROR : adding repeat region failed for $file"
		exit 1;
	fi    
	### intersect with miRbase bed file
	cat $file.sift.codons.map.repeat | awk 'NR>2' | awk '{print $1"\t"($2-1)"\t"$2}' > $file.sift.codons.map.repeat.ChrPos.bed
	if [ -s $file.sift.codons.map.repeat.ChrPos.bed ]
	then
		$bed/intersectBed -a $file.sift.codons.map.repeat.ChrPos.bed -b $miRbase -c | awk '{print $NF}' > $file.sift.codons.map.repeat.ChrPos.bed.i
	else
		touch $file.sift.codons.map.repeat.ChrPos.bed.i
	fi
	echo -e "\nmiRbase" >> $file.sift.codons.map.repeat.ChrPos.bed.i.tmp
	cat $file.sift.codons.map.repeat.ChrPos.bed.i >> $file.sift.codons.map.repeat.ChrPos.bed.i.tmp
	rm $file.sift.codons.map.repeat.ChrPos.bed.i $file.sift.codons.map.repeat.ChrPos.bed
	paste $file.sift.codons.map.repeat $file.sift.codons.map.repeat.ChrPos.bed.i.tmp > $file.sift.codons.map.repeat.base
	rm $file.sift.codons.map.repeat.ChrPos.bed.i.tmp
	num_a=`cat $file.sift.codons.map.repeat.base | wc -l `
	if [ $num_a == $num ]
	then
	    rm $file.sift.codons.map.repeat
	else
	    echo "ERROR : adding miRbase failed for $file"
		exit 1;
	fi    
	### intersect with SSR SCS
	## SSR=SNP Suspect Reason
	## SCS=SNP Clinical Significance
	cat $file.sift.codons.map.repeat.base | awk 'NR>2' | awk '{print $1"\t"($2-1)"\t"$2}' > $file.sift.codons.map.repeat.base.ChrPos.bed
	cat $ssr | awk -v num=chr$which_chr '$1 == num' > $file.sift.codons.map.repeat.base.ChrPos.bed.ssr.i
	cat $scs | awk -v num=chr$which_chr '$1 == num' > $file.sift.codons.map.repeat.base.ChrPos.bed.scs.i
	cat $sao | awk -v num=chr$which_chr '$1 == num' > $file.sift.codons.map.repeat.base.ChrPos.bed.sao.i
	cat $build | awk -v num=chr$which_chr '$1 == num' > $file.sift.codons.map.repeat.base.ChrPos.bed.build.i
	for type in ssr scs sao build
	do
		perl $script_path/match.pl -b $file.sift.codons.map.repeat.base.ChrPos.bed -i $file.sift.codons.map.repeat.base.ChrPos.bed.$type.i -o $file.sift.codons.map.repeat.base.ChrPos.bed.$type -t $type
	done
	
	echo -e "\nSNP_SuspectRegion" > $file.sift.codons.map.repeat.base.ChrPos.bed.ssr.tmp
	echo -e "\nSNP_ClinicalSig" > $file.sift.codons.map.repeat.base.ChrPos.bed.scs.tmp
	echo -e "\nVariant_AlleleOrigin" > $file.sift.codons.map.repeat.base.ChrPos.bed.sao.tmp
	echo -e "\nFirst_dbSNP_Build" > $file.sift.codons.map.repeat.base.ChrPos.bed.build.tmp
	
	for type in ssr scs sao build
	do
		cat $file.sift.codons.map.repeat.base.ChrPos.bed.$type >> $file.sift.codons.map.repeat.base.ChrPos.bed.$type.tmp
	done
	
	paste $file.sift.codons.map.repeat.base $file.sift.codons.map.repeat.base.ChrPos.bed.ssr.tmp $file.sift.codons.map.repeat.base.ChrPos.bed.scs.tmp $file.sift.codons.map.repeat.base.ChrPos.bed.sao.tmp $file.sift.codons.map.repeat.base.ChrPos.bed.build.tmp > $file.sift.codons.map.repeat.base.snp
	rm $file.sift.codons.map.repeat.base.ChrPos.bed.ssr.tmp $file.sift.codons.map.repeat.base.ChrPos.bed.scs.tmp $file.sift.codons.map.repeat.base.ChrPos.bed.sao.tmp $file.sift.codons.map.repeat.base.ChrPos.bed.build.tmp
	rm $file.sift.codons.map.repeat.base.ChrPos.bed.ssr $file.sift.codons.map.repeat.base.ChrPos.bed.scs $file.sift.codons.map.repeat.base.ChrPos.bed.sao $file.sift.codons.map.repeat.base.ChrPos.bed.build
	rm $file.sift.codons.map.repeat.base.ChrPos.bed
	rm $file.sift.codons.map.repeat.base.ChrPos.bed.ssr.i $file.sift.codons.map.repeat.base.ChrPos.bed.scs.i $file.sift.codons.map.repeat.base.ChrPos.bed.sao.i $file.sift.codons.map.repeat.base.ChrPos.bed.build.i 
	num=`cat $file.sift.codons.map.repeat.base.snp | wc -l`
	if [ $num == $num_a ]
	then
	    rm $file.sift.codons.map.repeat.base
	else
	    echo "ERROR : adding SSR and SCS failed for $file"
		exit 1;
	fi    
        ### add ucsc tracks
	cat $file.sift.codons.map.repeat.base.snp | awk 'NR>2' | awk '{print $1"\t"($2-1)"\t"$2}' > $file.sift.codons.map.repeat.base.snp.ChrPos.bed
	for i in conservation regulation tfbs tss enhancer
	do
	    if [ -s $file.sift.codons.map.repeat.base.snp.ChrPos.bed ]
		then
			$bed/intersectBed -b $file.sift.codons.map.repeat.base.snp.ChrPos.bed -a $UCSC/${i}.bed -wb > $file.sift.codons.map.repeat.base.snp.ChrPos.bed.${i}
	    else
			touch $file.sift.codons.map.repeat.base.snp.ChrPos.bed.${i}
		fi
		perl $script_path/matching_ucsc_tracks.pl $file.sift.codons.map.repeat.base.snp.ChrPos.bed $file.sift.codons.map.repeat.base.snp.ChrPos.bed.${i} $file.sift.codons.map.repeat.base.snp.ChrPos.bed.${i}.txt ${i} $chr
	done
	paste $file.sift.codons.map.repeat.base.snp $file.sift.codons.map.repeat.base.snp.ChrPos.bed.conservation.txt $file.sift.codons.map.repeat.base.snp.ChrPos.bed.regulation.txt $file.sift.codons.map.repeat.base.snp.ChrPos.bed.tfbs.txt $file.sift.codons.map.repeat.base.snp.ChrPos.bed.tss.txt $file.sift.codons.map.repeat.base.snp.ChrPos.bed.enhancer.txt > $file.sift.codons.map.repeat.base.snp.UCSCtracks 
	rm $file.sift.codons.map.repeat.base.snp.ChrPos.bed
	for i in conservation regulation tfbs tss enhancer
	do
	    rm $file.sift.codons.map.repeat.base.snp.ChrPos.bed.${i}.txt
	    rm $file.sift.codons.map.repeat.base.snp.ChrPos.bed.${i}
	done
	num_a=`cat $file.sift.codons.map.repeat.base.snp.UCSCtracks | wc -l`
	if [ $num_a == $num ]
	then
	    rm $file.sift.codons.map.repeat.base.snp
	else
	    echo "ERROR : addding UCSC tracks failed for $file"
		exit 1;
	fi
	
	### add polyphen
	if [ -f $poly/$sample.variants.chr${which_chr}.SNV.filter.i.c.vcf.poly ]
	then
		perl $script_path/add_polphen.pl $file.sift.codons.map.repeat.base.snp.UCSCtracks $poly/$sample.variants.chr${which_chr}.SNV.filter.i.c.vcf.poly $file.sift.codons.map.repeat.base.snp.UCSCtracks.poly
	else
		$script_path/email.sh $poly/$sample.variants.chr${which_chr}.SNV.filter.i.c.vcf.poly  "not found" $JOB_NAME $JOB_ID $run_info
		touch $poly/$sample.variants.chr${which_chr}.SNV.filter.i.c.vcf.poly.fix.log
		while [ -f $poly/$sample.variants.chr${which_chr}.SNV.filter.i.c.vcf.poly.fix.log ]
		do
			echo "waiting for the fix"
			sleep 2m
		done
		perl $script_path/add_polphen.pl $file.sift.codons.map.repeat.base.snp.UCSCtracks $poly/$sample.variants.chr${which_chr}.SNV.filter.i.c.vcf.poly $file.sift.codons.map.repeat.base.snp.UCSCtracks.poly
	fi
	num=`cat $file.sift.codons.map.repeat.base.snp.UCSCtracks.poly | wc -l `
	if [ $num == $num_a ]
	then
	    rm $file.sift.codons.map.repeat.base.snp.UCSCtracks
	else
	    echo "ERROR : adding polyphen failed for $file"
	fi    
   
	### add snpeff prediction
	perl $script_path/add_snpeff.pl -i $file.sift.codons.map.repeat.base.snp.UCSCtracks.poly -s $snpeff/$sample.chr${which_chr}.snv.eff -o $TempReports/$sample.chr${which_chr}.SNV.report
	perl $script_path/add_snpeff.pl -i $file.sift.codons.map.repeat.base.snp.UCSCtracks.poly -s $snpeff/$sample.chr${which_chr}.snv.filtered.eff -o $TempReports/$sample.chr${which_chr}.filtered.SNV.report
	num_a=`cat $TempReports/$sample.chr${which_chr}.SNV.report | awk -F'\t' '{print $1"_"$2}' | sort | uniq | wc -l`
	num_b=`cat $TempReports/$sample.chr${which_chr}.filtered.SNV.report | wc -l `
    rm $file.sift.codons.map.repeat.base.snp.UCSCtracks.poly	
	### add the entrez id to the report
	for report in $TempReports/$sample.chr${which_chr}.SNV.report $TempReports/$sample.chr${which_chr}.filtered.SNV.report
	do
		perl $script_path/add_entrezID.pl -i $report -m $GeneIdMap -o $report.entrezid
		mv $report.entrezid $report
		## exclude the redundant columns from the report
		perl $script_path/to.exclude.redundant.columns.from.report.pl $report $report.formatted
		mv $report.formatted $report
	done
        
	### add igv, tissue, pathway, kaviar column to the reports
	perl $script_path/add.cols.pl $TempReports/$sample.chr${which_chr}.SNV.report $run_info SNV > $TempReports/$sample.chr${which_chr}.SNV.xls
	perl $script_path/add.cols.pl $TempReports/$sample.chr${which_chr}.filtered.SNV.report $run_info SNV > $TempReports/$sample.chr${which_chr}.filtered.SNV.xls
	rm $TempReports/$sample.chr${which_chr}.SNV.report
	rm $TempReports/$sample.chr${which_chr}.filtered.SNV.report
	echo `date`	
fi