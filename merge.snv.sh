#!/bin/sh

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
		
if [ $# != 7 ];
then
	echo "Usage:<Tempreports><sample name><chromosome><sseq dir><sift dir><variant file><run info>";
else	
	set -x
	echo `date`
	TempReports=$1
	sample=$2 
	which_chr=$3
	sift=$4 
	sseq=$5
	var=$6	
	run_info=$7
	tool_info=$( cat $run_info | grep -w '^TOOL_INFO' | cut -d '=' -f2)
	script_path=$( cat $tool_info | grep -w '^WHOLEGENOME_PATH' | cut -d '=' -f2 )
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
	typeset -i codon
	typeset -i SNP_Type
	
	### add sift columns
        file=$TempReports/$var.rsIDs.allele_frequencies
	chr=`awk -F '\t' '{ for(i=1;i<=NF;i++){ if ($i == "Chr") {print i} } }' $file`
	pos=`awk -F '\t' '{ for(i=1;i<=NF;i++){ if ($i == "Position") {print i} } }' $file`
	ref=`awk -F '\t' '{ for(i=1;i<=NF;i++){ if ($i == "Ref") {print i} } }' $file`
	alt=`awk -F '\t' '{ for(i=1;i<=NF;i++){ if ($i == "Alt") {print i} } }' $file`
	perl $script_path/parse_siftPredictions.pl -i $file -s $sift/${sample}_chr${which_chr}_predictions.tsv -c $chr -p $pos -r $ref -a $alt -o $file.sift
	
	### sort the erport using the start and position
	perl $script_path/sort.variantReport.pl -i $file.sift -o $file.sift.sort -f Position
	mv $file.sift.sort $file.sift
	##  add codon preference
	touch $file.sift.forCodons
	cat $file.sift >> $file.sift.forCodons
	codon=`awk -F '\t' '{ for(i=1;i<=NF;i++){ if ($i ~ "Codons") {print i} } }' $file.sift.forCodons`
	SNP_Type=`awk -F '\t' '{ for(i=1;i<=NF;i++){ if ($i ~ "SNP Type") {print i} } }' $file.sift.forCodons`
	cat $file.sift.forCodons | cut -f ${codon},${SNP_Type} > $file.sift.forCodons.2columns
	perl $script_path/codon_pref.pl $codon_ref $file.sift.forCodons.2columns $file.sift.forCodons.2columns.added
	dos2unix $file.sift.forCodons $file.sift.forCodons.2columns.added
	paste $file.sift.forCodons $file.sift.forCodons.2columns.added > $file.sift.codons
	rm $file.sift.forCodons
	rm $file.sift.forCodons.2columns.added
	rm $file.sift.forCodons.2columns
	
        ### blacklisted, Alignablility or Uniqueness columns
	cat $file.sift.codons | cut -f 1,2 > $file.sift.codons.ChrPos
	sed -i '1d;2d' $file.sift.codons.ChrPos
	cat $file.sift.codons.ChrPos | awk '{print $1"\t"($2-1)"\t"$2}' > $file.sift.codons.ChrPos.bed
	rm $file.sift.codons.ChrPos
	$bed/intersectBed -a $file.sift.codons.ChrPos.bed -b $blacklisted -c | awk '{print $NF}' > $file.sift.codons.ChrPos.bed.black.i
	echo -e "\nBlacklistedRegion" > $file.sift.codons.ChrPos.bed.black.i.tmp
	cat  $file.sift.codons.ChrPos.bed.black.i.tmp $file.sift.codons.ChrPos.bed.black.i > $file.sift.codons.ChrPos.bed.black
	rm $file.sift.codons.ChrPos.bed.black.i.tmp $file.sift.codons.ChrPos.bed.black.i
	$bed/intersectBed -a $file.sift.codons.ChrPos.bed -b $mapability -c | awk '{print $NF}' > $file.sift.codons.ChrPos.bed.map.i
	echo -e "\nAlignability/Uniquness" > $file.sift.codons.ChrPos.bed.map.i.tmp
	cat  $file.sift.codons.ChrPos.bed.map.i.tmp $file.sift.codons.ChrPos.bed.map.i > $file.sift.codons.ChrPos.bed.map
	rm $file.sift.codons.ChrPos.bed.map.i.tmp $file.sift.codons.ChrPos.bed.map.i
	paste $file.sift.codons $file.sift.codons.ChrPos.bed.black $file.sift.codons.ChrPos.bed.map > $file.sift.codons.map
	rm $file.sift.codons.ChrPos.bed $file.sift.codons.ChrPos.bed.black $file.sift.codons.ChrPos.bed.map
	
	## intersect with repeat region bed file
	cat $file.sift.codons.map | cut -f 1,2 > $file.sift.codons.map.ChrPos
	sed -i '1d;2d' $file.sift.codons.map.ChrPos
	cat $file.sift.codons.map.ChrPos | awk '{print $1"\t"($2-1)"\t"$2}' > $file.sift.codons.map.ChrPos.bed
	rm $file.sift.codons.map.ChrPos
	$bed/intersectBed -a $file.sift.codons.map.ChrPos.bed -b $repeat -c | awk '{print $NF}' > $file.sift.codons.map.ChrPos.bed.i
	echo -e "\nRepeat_Region" >> $file.sift.codons.map.ChrPos.bed.i.tmp
	cat $file.sift.codons.map.ChrPos.bed.i >> $file.sift.codons.map.ChrPos.bed.i.tmp
	paste $file.sift.codons.map $file.sift.codons.map.ChrPos.bed.i.tmp > $file.sift.codons.map.repeat
	rm $file.sift.codons.map.ChrPos.bed.i.tmp $file.sift.codons.map.ChrPos.bed.i $file.sift.codons.map.ChrPos.bed
	
	### intersect with miRbase bed file
	cat $file.sift.codons.map.repeat | cut -f 1,2 > $file.sift.codons.map.repeat.ChrPos
	sed -i '1d;2d' $file.sift.codons.map.repeat.ChrPos
	cat $file.sift.codons.map.repeat.ChrPos | awk '{print $1"\t"($2-1)"\t"$2}' > $file.sift.codons.map.repeat.ChrPos.bed
	rm $file.sift.codons.map.repeat.ChrPos
	$bed/intersectBed -a $file.sift.codons.map.repeat.ChrPos.bed -b $miRbase -c | awk '{print $NF}' > $file.sift.codons.map.repeat.ChrPos.bed.i
	echo -e "\nmiRbase" >> $file.sift.codons.map.repeat.ChrPos.bed.i.tmp
	cat $file.sift.codons.map.repeat.ChrPos.bed.i >> $file.sift.codons.map.repeat.ChrPos.bed.i.tmp
	rm $file.sift.codons.map.repeat.ChrPos.bed.i $file.sift.codons.map.repeat.ChrPos.bed
	paste $file.sift.codons.map.repeat $file.sift.codons.map.repeat.ChrPos.bed.i.tmp > $file.sift.codons.map.repeat.base
	rm $file.sift.codons.map.repeat.ChrPos.bed.i.tmp
	
	### intersect with SSR SCS
	## SSR=SNP Suspect Reason
	## SCS=SNP Clinical Significance
	cat $file.sift.codons.map.repeat.base | cut -f 1,2 > $file.sift.codons.map.repeat.base.ChrPos
	sed -i '1d;2d' $file.sift.codons.map.repeat.base.ChrPos
	cat $file.sift.codons.map.repeat.base.ChrPos | awk '{print $1"\t"($2-1)"\t"$2}' > $file.sift.codons.map.repeat.base.ChrPos.bed
	rm $file.sift.codons.map.repeat.base.ChrPos
	$bed/intersectBed -a $file.sift.codons.map.repeat.base.ChrPos.bed -b $ssr -wb | awk '{print $(NF-3)"\t"$(NF-2)"\t"$(NF-1)"\t"$NF}' > $file.sift.codons.map.repeat.base.ChrPos.bed.ssr.i
	$bed/intersectBed -a $file.sift.codons.map.repeat.base.ChrPos.bed -b $scs -wb | awk '{print $(NF-3)"\t"$(NF-2)"\t"$(NF-1)"\t"$NF}' > $file.sift.codons.map.repeat.base.ChrPos.bed.scs.i
	
	perl $script_path/match.pl -b $file.sift.codons.map.repeat.base.ChrPos.bed -i $file.sift.codons.map.repeat.base.ChrPos.bed.ssr.i -o $file.sift.codons.map.repeat.base.ChrPos.bed.ssr
	perl $script_path/match.pl -b $file.sift.codons.map.repeat.base.ChrPos.bed -i $file.sift.codons.map.repeat.base.ChrPos.bed.scs.i -o $file.sift.codons.map.repeat.base.ChrPos.bed.scs
	echo -e "\nSNP_SuspectRegion" >> $file.sift.codons.map.repeat.base.ChrPos.bed.ssr.tmp
	echo -e "\nSNP_ClinicalSig" >> $file.sift.codons.map.repeat.base.ChrPos.bed.scs.tmp
	cat $file.sift.codons.map.repeat.base.ChrPos.bed.ssr >> $file.sift.codons.map.repeat.base.ChrPos.bed.ssr.tmp
	cat $file.sift.codons.map.repeat.base.ChrPos.bed.scs >> $file.sift.codons.map.repeat.base.ChrPos.bed.scs.tmp
	paste $file.sift.codons.map.repeat.base $file.sift.codons.map.repeat.base.ChrPos.bed.ssr.tmp $file.sift.codons.map.repeat.base.ChrPos.bed.scs.tmp > $file.sift.codons.map.repeat.base.snp
	rm $file.sift.codons.map.repeat.base.ChrPos.bed.ssr.tmp $file.sift.codons.map.repeat.base.ChrPos.bed.scs.tmp
	rm $file.sift.codons.map.repeat.base.ChrPos.bed.ssr $file.sift.codons.map.repeat.base.ChrPos.bed.scs
	rm $file.sift.codons.map.repeat.base.ChrPos.bed
	rm $file.sift.codons.map.repeat.base.ChrPos.bed.ssr.i $file.sift.codons.map.repeat.base.ChrPos.bed.scs.i
	
    ### add ucsc tracks
	cat $file.sift.codons.map.repeat.base.snp | cut -f 1,2 > $file.sift.codons.map.repeat.base.snp.ChrPos
	sed -i '1d;2d' $file.sift.codons.map.repeat.base.snp.ChrPos
	cat $file.sift.codons.map.repeat.base.snp.ChrPos | awk '{print $1"\t"($2-1)"\t"$2}' > $file.sift.codons.map.repeat.base.snp.ChrPos.bed
	rm $file.sift.codons.map.repeat.base.snp.ChrPos
	for i in conservation regulation tfbs tss enhancer
	do
		$bed/intersectBed -b $file.sift.codons.map.repeat.base.snp.ChrPos.bed -a $UCSC/${i}.bed -wb > $file.sift.codons.map.repeat.base.snp.ChrPos.bed.${i}
		perl $script_path/matching_ucsc_tracks.pl $file.sift.codons.map.repeat.base.snp.ChrPos.bed $file.sift.codons.map.repeat.base.snp.ChrPos.bed.${i} $file.sift.codons.map.repeat.base.snp.ChrPos.bed.${i}.txt ${i} $chr
		dos2unix $file.sift.codons.map.repeat.base.snp.repeat.base.snp.ChrPos.bed.${i}.txt 
	done
	cat $file.sift.codons.map.repeat.base.snp | sed 's/[ \t]*$//' > $file.sift.codons.map.repeat.base.snp.forUCSCtracks
	paste $file.sift.codons.map.repeat.base.snp.forUCSCtracks $file.sift.codons.map.repeat.base.snp.ChrPos.bed.conservation.txt $file.sift.codons.map.repeat.base.snp.ChrPos.bed.regulation.txt $file.sift.codons.map.repeat.base.snp.ChrPos.bed.tfbs.txt $file.sift.codons.map.repeat.base.snp.ChrPos.bed.tss.txt $file.sift.codons.map.repeat.base.snp.ChrPos.bed.enhancer.txt > $file.sift.codons.map.repeat.base.snp.UCSCtracks  
	rm $file.sift.codons.map.repeat.base.snp.forUCSCtracks 
        rm $file.sift.codons.map.repeat.base.snp.ChrPos.bed
	for i in conservation regulation tfbs tss enhancer
	do
		rm $file.sift.codons.map.repeat.base.snp.ChrPos.bed.${i}.txt
		rm $file.sift.codons.map.repeat.base.snp.ChrPos.bed.${i}
	done
	
	### add sseq annnoattion to the report
	perl $script_path/MergeListReport_SeattleSeq.pl -i $file.sift.codons.map.repeat.base.snp.UCSCtracks -s $sseq/$sample.chr${which_chr}.snv.sseq -c $chr -p $pos -r $ref -a $alt -o $TempReports/$sample.chr${which_chr}.SNV.report
	### add the entrez id to the report
	report=$TempReports/$sample.chr${which_chr}.SNV.report
	perl $script_path/add_entrezID.pl -i $report -m $GeneIdMap -o $report.entrezid
	mv $report.entrezid $report
	## exclude the redundant columns from the report
	perl $script_path/to.exclude.redundant.columns.from.report.pl $report $report.formatted
	mv $report.formatted $report
	echo `date`	
fi