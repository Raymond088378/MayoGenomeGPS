#!/bin/bash

if [ $# != 4 ]
then
    echo -e "script to parse a vcf file to tab delimited file using vctools vcf-query\nUsage: ./parse.vcf.sh <in_vcf_file> <out_file> <run_info><type of variant>";
else
    set -x
    echo `date`
    file=$1
    outfile=$2
    run_info=$3
    type=$4
    tool_info=$( cat $run_info | grep -w '^TOOL_INFO' | cut -d '=' -f2)
    tabix=$( cat $tool_info | grep -w '^TABIX' | cut -d '=' -f2 )
    vcftools=$( cat $tool_info | grep -w '^VCFTOOLS' | cut -d '=' -f2 )
	perllib=$( cat $tool_info | grep -w '^PERLLIB_VCF' | cut -d '=' -f2)
	script_path=$( cat $tool_info | grep -w '^WORKFLOW_PATH' | cut -d '=' -f2 )
	export PERL5LIB=$perllib:$PERL5LIB
	export PATH=$tabix/:$PATH
	file_name=`basename $file`	
	dir=`dirname $outfile`	
    num=`cat $file | grep -v '#' | wc -l`
	cat $file | $script_path/convert.vcf.pl > $dir/$file_name

    $tabix/bgzip $dir/$file_name
    $tabix/tabix -p vcf $dir/$file_name.gz

    
    sample_cols=""
	samples=""
    for i in `$vcftools/bin/vcf-query -l $dir/$file_name.gz` 
    do 
        if [ $type == 'SNV' ]
        then
            samples="$samples""$i"'::::::' 
            sample_cols=':GenotypeClass:Ref-SupportedReads:Alt-SupportedReads:ReadDepth:Quality:CloseToIndel'$sample_cols;
        else
            samples="$samples""$i"'::';
            sample_cols=':Indel-supportedRead:ReadDepth'$sample_cols;
        fi    
    done

    samples=`echo $samples | sed -e 's/\(.*\)./\1/' | tr ':' '\t'`;
    sample_cols=`echo $sample_cols | tr ':' '\t'`;

    if [ $type == 'SNV' ]
    then
        echo -e "\t\t\t\t\t\t$samples">$outfile
        echo -e "Chr\tPosition\tInCaptureKit\t#AlternateHits\tRef\tAlt${sample_cols}" >>$outfile
		$vcftools/bin/vcf-query $dir/$file_name.gz -f '%CHROM\t%POS\t%INFO/CAPTURE\t%INFO/ED\t%REF\t%ALT\t[%GT\t%AD\t%DP\t%GQ\t%C2I\t]\n'	| $script_path/parse.snv.vcftools.pl >> $outfile
    else
        echo -e "\t\t\t\t\t\t\t\t$samples">$outfile
        echo -e "Chr\tStart\tStop\tInCaptureKit\t#AlternateHits\tRef\tAlt\tBase-Length${sample_cols}" >>$outfile
        $vcftools/bin/vcf-query $dir/$file_name.gz -f '%CHROM\t%POS\t%INFO/CAPTURE\t%INFO/ED\t%REF\t%ALT\t[%AD\t%DP\t]\n' | $script_path/parse.pl >> $outfile
    fi    
    num_a=`cat $outfile | awk 'NR>2' | wc -l`
	if [ $num_a == $num ]
	then
		rm $dir/$file_name.gz.tbi $dir/$file_name.gz
    else
		touch $outfile.fix.log
		$script_path/email.sh $outfile "not found" OnTarget_variant.sh $run_info
		$script_path/wait.sh $outfile.fix.log
		rm $dir/$file_name.gz.tbi $dir/$file_name.gz
	fi
	echo `date`
fi