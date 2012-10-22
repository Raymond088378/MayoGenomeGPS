#!/bin/bash

if [ $# != 7 ]
then
    echo -e "Usage: stand alone script for annotation using a single vcf or text file \n <sample name><vcf or txt file><output folder><tool_info file> <genome build><single/multi><script_path>"
else
	echo -e "\n ************************************************** \n"
	echo " Started the annotation for your data on `date` "
	sample=$1
	file=$2
	output=$3
	tool_info=$4
	GenomeBuild=$5
	multi=$6
	script_path=$7
	## paths
	mkdir -p $output
	cp $tool_info $output/$sample.tool_info.txt
	tool_info=$output/$sample.tool_info.txt
	snpeff=$( cat $tool_info | grep -w '^SNPEFF' | cut -d '=' -f2)
	sift=$( cat $tool_info | grep -w '^SIFT' | cut -d '=' -f2)
	pph=$(cat $tool_info | grep -w '^POLYPHEN' |  cut -d '=' -f2)
	sift_ref=$( cat $tool_info | grep -w '^SIFT_REF' | cut -d '=' -f2)
	samtools=$( cat $tool_info | grep -w '^SAMTOOLS' | cut -d '=' -f2 )
	ref=$( cat $tool_info | grep -w '^REF_GENOME' | cut -d '=' -f2)
	gatk=$( cat $tool_info | grep -w '^GATK' | cut -d '=' -f2)
	vcftools=$( cat $tool_info | grep -w '^VCFTOOLS' | cut -d '=' -f2)
	perllib=$( cat $tool_info | grep -w '^PERLLIB_VCF' | cut -d '=' -f2)
	thread=$( cat $tool_info | grep -w '^THREADS' | cut -d '=' -f2)
    java=$( cat $tool_info | grep -w '^JAVA' | cut -d '=' -f2)
	http=$( cat $tool_info | grep -w '^HTTP_SERVER' | cut -d '=' -f2)
	export PERL5LIB=/projects/bsi/bictools/apps/annotation/polyphen/2.2.2/perl/:$perllib:$PERL5LIB
	tabix=$( cat $tool_info | grep -w '^TABIX' | cut -d '=' -f2 )
    vcftools=$( cat $tool_info | grep -w '^VCFTOOLS' | cut -d '=' -f2 )
	perllib=$( cat $tool_info | grep -w '^PERLLIB_VCF' | cut -d '=' -f2)
	blat_params=$( cat $tool_info | grep -w '^BLAT_params' | cut -d '=' -f2 )
	export PERL5LIB=$perllib:$PERL5LIB
	PATH=$tabix/:$PATH
	echo " vcf validation step "
	type=`cat $file | head -1 | awk '{if ($0 ~ /^##/) print "vcf";else print "txt"}'`
	ff="$sample.vcf"
	if [ $type == "txt" ]
	then
		if [[ `cat $file| awk -F'\t' '{print NF}' | sort | uniq` != 4 ]]
		then
			echo "txt file is not properly formatted"
			exit 1;
		else	
			$script_path/convert_txt_vcf.pl $file $sample > $output/$ff
		fi
	else
		if [ $multi == "multi" ]
		then
			cat $file  | $script_path/convert.vcf.pl > $output/$ff
		else
			col=`cat $file | awk '$0 ~ /#CHROM/' | awk -v num=$sample -F '\t' '{ for(i=1;i<=NF;i++){ if ($i == num) {print i} } }'`
			cat $file | awk -v num=$col 'BEGIN {OFS="\t"} { if ($0 ~ /^##/) print $0; else print $1,$2,$3,$4,$5,$6,$7,$8,$9,$num; }' | $script_path/convert.vcf.pl > $output/$ff
		fi
	fi
	
	### extract multi-allelic variants and keep that as a vcf file
	cat $output/$ff | awk '$0 ~ /^#/ || $5 ~ /,/ || $4 ~ /,/' > $output/$sample.multi.vcf
	cat $output/$ff | awk '$0 ~ /^#/ || ($5 !~ /,/ && $4 !~ /,/)' > $output/$ff.temp
	mv $output/$ff.temp $output/$ff
	
	## ADD BLAT column
	n=`cat $output/$ff |  awk '$0 ~ /^#/' | awk '$0 ~ /^##INFO=<ID=ED/' | wc -l`
	if [ $n == 0 ]
	then
		echo " Adding Blat column to the vcf file"
        blat=$( cat $tool_info | grep -w '^BLAT' | cut -d '=' -f2 )
		blat_ref=$( cat $tool_info | grep -w '^BLAT_REF' | cut -d '=' -f2 )
		$script_path/vcf_blat_verify.pl -i $output/$ff -o $output/$ff.tmp -b $blat -r $ref -sam $samtools -br $blat_ref $blat_params
		mv $output/$ff.tmp $output/$ff
	fi
	$script_path/vcfsort.pl $ref.fai $output/$ff > $output/$ff.sort
	mv $output/$ff.sort $output/$ff
	n=`cat $output/$ff |   awk '$0 ~ /^#/' | awk '$0 ~ /^##INFO=<ID=CAPTURE/' | wc -l`
	perl $script_path/vcf_to_variant_vcf.pl -i $output/$ff -v $output/$ff.SNV.vcf -l $output/$ff.INDEL.vcf -t both
	rm $output/$ff

	if [ $n == 0 ]
	then
		echo " Adding CloseToIndel column to the vcf file"
		$script_path/markSnv_IndelnPos.pl -s $output/$ff.SNV.vcf -i $output/$ff.INDEL.vcf -n 10 -o $output/$ff.SNV.vcf.pos
		cat $output/$ff.SNV.vcf.pos | $script_path/add.info.close2indel.vcf.pl | $script_path/add.info.capture.vcf.pl > $output/$ff.SNV.vcf
		rm $output/$ff.SNV.vcf.pos
	fi
	## RUN SIFT
	echo " Running SIFT "
	$script_path/sift.sh $output $ff $sift $sift_ref $script_path $sample $thread &

	### RUN SNPEFF
	echo " Running SNPEFF "
	$script_path/snpeff.sh $java $snpeff $GenomeBuild $output $ff $gatk $ref $vcftools $script_path $sample &
	### POLYPHEN
	echo " Running Polyphen "
	$script_path/polyphen.sh $output $ff $pph $script_path $GenomeBuild $sample $thread &


	#### now add annotations
	bed=$( cat $tool_info | grep -w '^BEDTOOLS' | cut -d '=' -f2 )
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
	dbsnp_rsids_snv=$( cat $tool_info | grep -w '^dbSNP_SNV_rsIDs' | cut -d '=' -f2)
	dbsnp_rsids_disease=$( cat $tool_info | grep -w '^dbSNP_disease_rsIDs' | cut -d '=' -f2)
	hapmap=$( cat $tool_info | grep -w '^HAPMAP' | cut -d '=' -f2)
	kgenome=$( cat $tool_info | grep -w '^KGENOME' | cut -d '=' -f2)
	bgi=$( cat $tool_info | grep -w '^BGI_REF' | cut -d '=' -f2)
	cosmic=$( cat $tool_info | grep -w '^COSMIC_SNV_REF' | cut -d '=' -f2)
	esp=$( cat $tool_info | grep -w '^ESP' | cut -d '=' -f2)

	typeset -i codon
	typeset -i SNP_Type
	#convert to text file
	echo " Converting VCF files to TEXT files"
    $script_path/parse.vcf.sh $output/$ff.SNV.vcf $output/$sample.snv $tool_info SNV $script_path $sample
	$script_path/parse.vcf.sh $output/$ff.INDEL.vcf $output/$sample.indel $tool_info INDEL $script_path $sample

	file=$output/$sample.snv
	## Add DBSNP
	echo " Adding annotation columns to excel worksheet "
	cat $file | awk 'NR>1' > $file.forrsIDs
	$script_path/add.rsids.pl -i $file.forrsIDs -s $dbsnp_rsids_snv -o $file.forrsIDs.added
	$script_path/add.dbsnp.disease.snv.pl -i $file.forrsIDs.added -b 1 -s $dbsnp_rsids_disease -c 1 -p 2 -o $file.forrsIDs.added.disease
	$script_path/extract.rsids.pl -i $file -r $file.forrsIDs.added.disease -o $file.rsIDs -v SNV
	rm $file.forrsIDs.added $file.forrsIDs.added.disease $file.forrsIDs $file
	echo " dbSNP columns added to the report "

	### Add Frequencies
	cat $file.rsIDs | awk 'NR>1' | cut -f 1-5 > $file.rsIDs.forFrequencies
	echo " adding allele frequencies to the report "
	$script_path/add_hapmap_1kgenome_allele_freq.pl -i $file.rsIDs.forFrequencies -c 1 -p 2 -b 1 -e CEU -s $hapmap/all_allele_freqs_CEU.txt -g $kgenome/CEU.$GenomeBuild -o $file.rsIDs.CEU&&perl $script_path/add_hapmap_1kgenome_allele_freq.pl -i $file.rsIDs.CEU -c 1 -p 2 -b 1 -e YRI -s $hapmap/all_allele_freqs_YRI.txt -g $kgenome/YRI.$GenomeBuild -o $file.rsIDs.CEU.YRI&&perl $script_path/add_hapmap_1kgenome_allele_freq.pl -i $file.rsIDs.CEU.YRI -c 1 -p 2 -b 1 -e JPT+CHB -s $hapmap/all_allele_freqs_JPT+CHB.txt -g $kgenome/JPT+CHB.$GenomeBuild -o $file.rsIDs.CEU.YRI.CHBJPT.txt
	rm $file.rsIDs.CEU $file.rsIDs.CEU.YRI
	## BGI
	perl $script_path/add_bgi_freq.pl -i $file.rsIDs.CEU.YRI.CHBJPT.txt -r $bgi -o $file.rsIDs.CEU.YRI.CHBJPT.BGI.txt
	rm $file.rsIDs.CEU.YRI.CHBJPT.txt
	## Cosmic
	perl $script_path/add.cosmic.pl $file.rsIDs.CEU.YRI.CHBJPT.BGI.txt 1 $cosmic $GenomeBuild 1 $file.rsIDs.CEU.YRI.CHBJPT.BGI.Cosmic.txt
	rm $file.rsIDs.CEU.YRI.CHBJPT.BGI.txt
	## ESP
	perl $script_path/add_esp.pl -i $file.rsIDs.CEU.YRI.CHBJPT.BGI.Cosmic.txt -d $esp -o $file.rsIDs.CEU.YRI.CHBJPT.BGI.Cosmic.ESP.txt
	perl $script_path/extract.allele_freq.pl -i $file.rsIDs -f $file.rsIDs.CEU.YRI.CHBJPT.BGI.Cosmic.ESP.txt -o $file.rsIDs.allele_frequencies -v SNV
	rm $file.rsIDs.CEU.YRI.CHBJPT.BGI.Cosmic.txt $file.rsIDs.CEU.YRI.CHBJPT.BGI.Cosmic.ESP.txt $file.rsIDs.forFrequencies $file.rsIDs
	echo " Allele frequencies added to the report "

	file=$file.rsIDs.allele_frequencies
	chr=`awk -F '\t' '{ for(i=1;i<=NF;i++){ if ($i == "Chr") {print i} } }' $file`
	pos=`awk -F '\t' '{ for(i=1;i<=NF;i++){ if ($i == "Position") {print i} } }' $file`
	ref=`awk -F '\t' '{ for(i=1;i<=NF;i++){ if ($i == "Ref") {print i} } }' $file`
	alt=`awk -F '\t' '{ for(i=1;i<=NF;i++){ if ($i == "Alt") {print i} } }' $file`

	while [[ ! -s $output/$sample.predictions.tsv ]]
	do
		echo "waiting for sift to finish"
		sleep 2m
	done

	## sift
	perl $script_path/parse_siftPredictions.pl -i $file -s $output/$sample.predictions.tsv -c $chr -p $pos -r $ref -a $alt -o $file.sift
	rm $output/$sample.predictions.tsv $file
	perl $script_path/sort.variantReport.pl -i $file.sift -o $file.sift.sort -f Position
	mv $file.sift.sort $file.sift
	echo " sorted report is ready to add more annotation "
	##  add codon preference
	codon=`awk -F '\t' '{ for(i=1;i<=NF;i++){ if ($i ~ "Codons") {print i} } }' $file.sift`
	SNP_Type=`awk -F '\t' '{ for(i=1;i<=NF;i++){ if ($i ~ "SNP Type") {print i} } }' $file.sift`
	cat $file.sift | cut -f ${codon},${SNP_Type} > $file.sift.forCodons.2columns
	perl $script_path/codon_pref.pl $codon_ref $file.sift.forCodons.2columns $file.sift.forCodons.2columns.added
	paste $file.sift $file.sift.forCodons.2columns.added > $file.sift.codons
	rm $file.sift.forCodons.2columns.added $file.sift.forCodons.2columns $file.sift
	echo " codon preference column is added "
	### blacklisted, Alignablility or Uniqueness columns
	cat $file.sift.codons | awk 'NR>2' | awk '{print $1"\t"($2-1)"\t"$2}' > $file.sift.codons.ChrPos.bed
	$bed/intersectBed -a $file.sift.codons.ChrPos.bed -b $blacklisted -c | awk '{print $NF}' > $file.sift.codons.ChrPos.bed.black.i
	echo -e "\nBlacklistedRegion" > $file.sift.codons.ChrPos.bed.black.i.tmp
	cat  $file.sift.codons.ChrPos.bed.black.i.tmp $file.sift.codons.ChrPos.bed.black.i > $file.sift.codons.ChrPos.bed.black
	rm $file.sift.codons.ChrPos.bed.black.i.tmp $file.sift.codons.ChrPos.bed.black.i
	$bed/intersectBed -a $file.sift.codons.ChrPos.bed -b $mapability -c | awk '{print $NF}' > $file.sift.codons.ChrPos.bed.map.i
	echo -e "\nAlignability/Uniquness" > $file.sift.codons.ChrPos.bed.map.i.tmp
	cat  $file.sift.codons.ChrPos.bed.map.i.tmp $file.sift.codons.ChrPos.bed.map.i > $file.sift.codons.ChrPos.bed.map
	rm $file.sift.codons.ChrPos.bed.map.i.tmp $file.sift.codons.ChrPos.bed.map.i
	paste $file.sift.codons $file.sift.codons.ChrPos.bed.black $file.sift.codons.ChrPos.bed.map > $file.sift.codons.map
	rm $file.sift.codons.ChrPos.bed $file.sift.codons.ChrPos.bed.black $file.sift.codons.ChrPos.bed.map $file.sift.codons
	## intersect with repeat region bed file
	echo " adding filtering tracks to the report "
	cat $file.sift.codons.map | awk 'NR>2' | awk '{print $1"\t"($2-1)"\t"$2}' > $file.sift.codons.map.ChrPos.bed
	$bed/intersectBed -a $file.sift.codons.map.ChrPos.bed -b $repeat -c | awk '{print $NF}' > $file.sift.codons.map.ChrPos.bed.i
	echo -e "\nRepeat_Region" >> $file.sift.codons.map.ChrPos.bed.i.tmp
	cat $file.sift.codons.map.ChrPos.bed.i >> $file.sift.codons.map.ChrPos.bed.i.tmp
	paste $file.sift.codons.map $file.sift.codons.map.ChrPos.bed.i.tmp > $file.sift.codons.map.repeat
	rm $file.sift.codons.map.ChrPos.bed.i.tmp $file.sift.codons.map.ChrPos.bed.i $file.sift.codons.map.ChrPos.bed $file.sift.codons.map
	### intersect with miRbase bed file
	echo " adding more columns "
        cat $file.sift.codons.map.repeat | awk 'NR>2' | awk '{print $1"\t"($2-1)"\t"$2}' > $file.sift.codons.map.repeat.ChrPos.bed
	$bed/intersectBed -a $file.sift.codons.map.repeat.ChrPos.bed -b $miRbase -c | awk '{print $NF}' > $file.sift.codons.map.repeat.ChrPos.bed.i
	echo -e "\nmiRbase" >> $file.sift.codons.map.repeat.ChrPos.bed.i.tmp
	cat $file.sift.codons.map.repeat.ChrPos.bed.i >> $file.sift.codons.map.repeat.ChrPos.bed.i.tmp
	rm $file.sift.codons.map.repeat.ChrPos.bed.i $file.sift.codons.map.repeat.ChrPos.bed
	paste $file.sift.codons.map.repeat $file.sift.codons.map.repeat.ChrPos.bed.i.tmp > $file.sift.codons.map.repeat.base
	rm $file.sift.codons.map.repeat.ChrPos.bed.i.tmp $file.sift.codons.map.repeat
	cat $file.sift.codons.map.repeat.base | awk 'NR>2' | awk '{print $1"\t"($2-1)"\t"$2}' > $file.sift.codons.map.repeat.base.ChrPos.bed
	ln -s $ssr $file.sift.codons.map.repeat.base.ChrPos.bed.ssr.i
	ln -s $scs $file.sift.codons.map.repeat.base.ChrPos.bed.scs.i
	ln -s $sao $file.sift.codons.map.repeat.base.ChrPos.bed.sao.i
	ln -s $build $file.sift.codons.map.repeat.base.ChrPos.bed.build.i
	echo " filtering tracks are added to the report "
	for type in ssr scs sao build
	do
		perl $script_path/match.pl -b $file.sift.codons.map.repeat.base.ChrPos.bed -i $file.sift.codons.map.repeat.base.ChrPos.bed.$type.i -o $file.sift.codons.map.repeat.base.ChrPos.bed.$type -t $type
	done
	echo -e "\nSNP_SuspectRegion" > $file.sift.codons.map.repeat.base.ChrPos.bed.ssr.tmp
	echo -e "\nSNP_ClinicalSig" > $file.sift.codons.map.repeat.base.ChrPos.bed.scs.tmp
	echo -e "\nVariant_AlleleOrigin" > $file.sift.codons.map.repeat.base.ChrPos.bed.sao.tmp
	echo -e "\nFirst_dbSNP_Build" > $file.sift.codons.map.repeat.base.ChrPos.bed.build.tmp
	echo " adding more columns "
        for type in ssr scs sao build
	do
		cat $file.sift.codons.map.repeat.base.ChrPos.bed.$type >> $file.sift.codons.map.repeat.base.ChrPos.bed.$type.tmp
	done

	paste $file.sift.codons.map.repeat.base $file.sift.codons.map.repeat.base.ChrPos.bed.ssr.tmp $file.sift.codons.map.repeat.base.ChrPos.bed.scs.tmp $file.sift.codons.map.repeat.base.ChrPos.bed.sao.tmp $file.sift.codons.map.repeat.base.ChrPos.bed.build.tmp > $file.sift.codons.map.repeat.base.snp
	rm $file.sift.codons.map.repeat.base.ChrPos.bed.ssr.tmp $file.sift.codons.map.repeat.base.ChrPos.bed.scs.tmp $file.sift.codons.map.repeat.base.ChrPos.bed.sao.tmp $file.sift.codons.map.repeat.base.ChrPos.bed.build.tmp
	rm $file.sift.codons.map.repeat.base.ChrPos.bed.ssr $file.sift.codons.map.repeat.base.ChrPos.bed.scs $file.sift.codons.map.repeat.base.ChrPos.bed.sao $file.sift.codons.map.repeat.base.ChrPos.bed.build
	rm $file.sift.codons.map.repeat.base.ChrPos.bed
	rm $file.sift.codons.map.repeat.base.ChrPos.bed.ssr.i $file.sift.codons.map.repeat.base.ChrPos.bed.scs.i $file.sift.codons.map.repeat.base.ChrPos.bed.sao.i $file.sift.codons.map.repeat.base.ChrPos.bed.build.i $file.sift.codons.map.repeat.base
	echo " adding more columns "
	cat $file.sift.codons.map.repeat.base.snp | awk 'NR>2' | awk '{print $1"\t"($2-1)"\t"$2}' > $file.sift.codons.map.repeat.base.snp.ChrPos.bed
	for i in conservation regulation tfbs tss enhancer
	do
		$bed/intersectBed -b $file.sift.codons.map.repeat.base.snp.ChrPos.bed -a $UCSC/${i}.bed -wb > $file.sift.codons.map.repeat.base.snp.ChrPos.bed.${i}
		perl $script_path/matching_ucsc_tracks.pl $file.sift.codons.map.repeat.base.snp.ChrPos.bed $file.sift.codons.map.repeat.base.snp.ChrPos.bed.${i} $file.sift.codons.map.repeat.base.snp.ChrPos.bed.${i}.txt ${i}
	done
	echo " UCSC tracks are added to the report "
	paste $file.sift.codons.map.repeat.base.snp $file.sift.codons.map.repeat.base.snp.ChrPos.bed.conservation.txt $file.sift.codons.map.repeat.base.snp.ChrPos.bed.regulation.txt $file.sift.codons.map.repeat.base.snp.ChrPos.bed.tfbs.txt $file.sift.codons.map.repeat.base.snp.ChrPos.bed.tss.txt $file.sift.codons.map.repeat.base.snp.ChrPos.bed.enhancer.txt > $file.sift.codons.map.repeat.base.snp.UCSCtracks
	rm $file.sift.codons.map.repeat.base.snp.ChrPos.bed $file.sift.codons.map.repeat.base.snp
	for i in conservation regulation tfbs tss enhancer
	do
	    rm $file.sift.codons.map.repeat.base.snp.ChrPos.bed.${i}.txt
	    rm $file.sift.codons.map.repeat.base.snp.ChrPos.bed.${i}
	done
	while [[ ! -f $output/$sample.polyphen.txt ||  ! -s $output/$sample.INDEL.filtered.eff || ! -s $output/$sample.SNV.filtered.eff ]]
	do
		echo " waiting for snpeff, polyphen annotation to complete "
		sleep 2m
	done
	rm $output/$ff.SNV.vcf $output/$ff.INDEL.vcf
	rm $output/log
	## polyphen
	perl $script_path/add_polphen.pl $file.sift.codons.map.repeat.base.snp.UCSCtracks $output/$sample.polyphen.txt $file.sift.codons.map.repeat.base.snp.UCSCtracks.poly
	rm $output/$sample.polyphen.txt $file.sift.codons.map.repeat.base.snp.UCSCtracks
	perl $script_path/add_snpeff.pl -i $file.sift.codons.map.repeat.base.snp.UCSCtracks.poly -s $output/$sample.SNV.eff -o $output/$sample.SNV.report
	perl $script_path/add_snpeff.pl -i $file.sift.codons.map.repeat.base.snp.UCSCtracks.poly -s $output/$sample.SNV.filtered.eff -o $output/$sample.filtered.SNV.report
	rm $output/$sample.SNV.eff $output/$sample.SNV.filtered.eff $file.sift.codons.map.repeat.base.snp.UCSCtracks.poly
	echo " polyphen and snpeff are added to the report "
        for report in $output/$sample.SNV.report $output/$sample.filtered.SNV.report
	do
		perl $script_path/add_entrezID.pl -i $report -m $GeneIdMap -o $report.entrezid
		mv $report.entrezid $report
		## exclude the redundant columns from the report
		perl $script_path/to.exclude.redundant.columns.from.report.pl $report $report.formatted
		mv $report.formatted $report
	done
	perl $script_path/add.cols.pl $output/$sample.SNV.report SNV $http $GenomeBuild > $output/$sample.SNV.xls
	perl $script_path/add.cols.pl $output/$sample.filtered.SNV.report SNV $http $GenomeBuild > $output/$sample.filtered.SNV.xls
	rm $output/$sample.SNV.report $output/$sample.filtered.SNV.report
        perl $script_path/sort.variantReport.pl -i $output/$sample.SNV.xls -o $output/$sample.SNV.xls.sort -f Position
        mv $output/$sample.SNV.xls.sort $output/$sample.SNV.xls
        perl $script_path/sort.variantReport.pl -i $output/$sample.filtered.SNV.xls -o $output/$sample.filtered.SNV.xls.sort -f Position
        mv $output/$sample.filtered.SNV.xls.sort $output/$sample.filtered.SNV.xls
	echo " SNV report is ready "
	##INDEL
	## DBSNP
	echo " adding indel annotation to the INDELs"
        dbsnp_rsids_indel=$( cat $tool_info | grep -w '^dbSNP_INDEL_rsIDs' | cut -d '=' -f2)
	file=$output/$sample.indel
	cat $file | awk 'NR>1' > $file.forrsIDs
	perl $script_path/add_dbsnp_indel.pl -i $file.forrsIDs -b 1 -s $dbsnp_rsids_indel -c 1 -p 2 -x 3 -o $file.forrsIDs.added
	cat $file.forrsIDs.added | awk '{if(NR != 1) print $0"\t0"; else print $0"\tDiseaseVariant"}' > $file.forrsIDs.added.disease
	perl $script_path/extract.rsids.pl -i $file -r $file.forrsIDs.added.disease -o $file.rsIDs -v INDEL
	rm $file.forrsIDs.added $file.forrsIDs.added.disease $file $file.forrsIDs
	## Frequency
        echo " adding more columns "
	cosmic=$( cat $tool_info | grep -w '^COSMIC_INDEL_REF' | cut -d '=' -f2)
	cat $file.rsIDs | awk 'NR>1' | cut -f 1-5 > $file.rsIDs.forfrequencies
	perl $script_path/add.cosmic.pl $file.rsIDs.forfrequencies 0 $cosmic $GenomeBuild 1 $file.rsIDs.forfrequencies.cosmic.txt
	rm $file.rsIDs.forfrequencies
	perl $script_path/extract.allele_freq.pl -i $file.rsIDs -f $file.rsIDs.forfrequencies.cosmic.txt -o $file.rsIDs.frequencies -v INDEL
	rm $file.rsIDs $file.rsIDs.forfrequencies.cosmic.txt
	## add snpeff prediction
	echo " annotated reports should be done soon "
	perl $script_path/add_snpeff_indel.pl -i $file.rsIDs.frequencies -s $output/$sample.INDEL.eff -o $output/$sample.INDEL.report
	perl $script_path/add_snpeff_indel_filter.pl -i $file.rsIDs.frequencies -s $output/$sample.INDEL.filtered.eff -o $output/$sample.filtered.INDEL.report
	rm $output/$sample.INDEL.eff $output/$sample.INDEL.filtered.eff $file.rsIDs.frequencies
	for report in $output/$sample.INDEL.report $output/$sample.filtered.INDEL.report
	do
		perl $script_path/add_entrezID.pl -i $report -m $GeneIdMap -o $report.entrezid
		mv $report.entrezid $report
		perl $script_path/to.exclude.redundant.columns.from.report.pl $report $report.formatted
		mv $report.formatted $report
	done
	perl $script_path/add.cols.pl $output/$sample.INDEL.report INDEL $http $GenomeBuild > $output/$sample.INDEL.xls
	perl $script_path/add.cols.pl $output/$sample.filtered.INDEL.report INDEL $http $GenomeBuild > $output/$sample.filtered.INDEL.xls
        rm $output/$sample.INDEL.report $output/$sample.filtered.INDEL.report
        perl $script_path/sort.variantReport.pl -i $output/$sample.INDEL.xls -o $output/$sample.INDEL.xls.sort -f Start
        mv $output/$sample.INDEL.xls.sort $output/$sample.INDEL.xls
        perl $script_path/sort.variantReport.pl -i $output/$sample.filtered.INDEL.xls -o $output/$sample.filtered.INDEL.xls.sort -f Start
        mv $output/$sample.filtered.INDEL.xls.sort $output/$sample.filtered.INDEL.xls
        echo " Annotation finished on `date` "
	echo -e "\n ************************************************** \n"
fi
