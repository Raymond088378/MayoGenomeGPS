output=$1
ff=$2
sift=$3
sift_ref=$4
script_path=$5
sample=$6
thread=$7

if [ `cat $output/$ff.SNV.vcf | awk '$0 !~ /^#/' | wc -l` -gt 0 ]
then
	cat $output/$ff.SNV.vcf | awk '$0 !~ /^#/' | sed -e '/chr/s///g' | awk '{print $1","$2",1,"$4"/"$5}' > $output/$ff.SNV.vcf.sift
	a=`pwd`

	#running SIFT for each sample
	cd $sift
	$script_path/parallel.sift.pl $thread $output/$ff.SNV.vcf.sift $sift_ref $sift $output/$sample.predictions.tsv $output/ $sample
	perl $script_path/sift.inconsistent.pl $output/$sample.predictions.tsv $output/$ff.SNV.vcf.sift
	mv $output/$sample.predictions.tsv_mod $output/$sample.predictions.tsv
	rm $output/$ff.SNV.vcf.sift
	cd $a
else
	 echo -e "Coordinates\tCodons\tTranscript ID\tProtein ID\tSubstitution\tRegion\tdbSNP ID\tSNP Type\tPrediction\tScore\tMedian Info\t# Seqs at position\tGene ID\tGene Name\tOMIM Disease\tAverage Allele Freqs\tUser Comment" >$output/$sample.predictions.tsv
fi	 
	echo " SIFT is done "
