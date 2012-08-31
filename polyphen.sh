#!/bin/bash

if [ $# -le 3 ]
then
    echo "Usage: wrapper to run polyphen \n <polphen dir> <input directory><sample><run info>"
else
    set -x
    echo `date`
    polyphen=$1
    input=$2
    sample=$3
    run_info=$4
    if [ $5 ]
    then
        prefix=$5
    fi	
	#SGE_TASK_ID=22
    chr=$(cat $run_info | grep -w '^CHRINDEX' | cut -d '=' -f2 | tr ":" "\n" | head -n $SGE_TASK_ID | tail -n 1)
    tool_info=$(cat $run_info | grep -w '^TOOL_INFO' |  cut -d '=' -f2)
    pph=$(cat $tool_info | grep -w '^POLYPHEN' |  cut -d '=' -f2)
    perl_lib=$(cat $tool_info | grep -w '^PERL_POLYPHEN_LIB' |  cut -d '=' -f2)
    genome_build=$(cat $run_info | grep -w '^GENOMEBUILD' |  cut -d '=' -f2)
    perl=$(cat $tool_info | grep -w '^PERL_CIRCOS' |  cut -d '=' -f2)
    script_path=$( cat $tool_info | grep -w '^WHOLEGENOME_PATH' | cut -d '=' -f2 )	
	threads=2   
	if [ $5 ]
    then
		sam=$prefix.$sample
    else
		sam=$sample
    fi
    
    export PERL5LIB=$perl_lib:${PERL5LIB}
    
    snv_file=$sam.variants.chr$chr.SNV.filter.i.c.vcf
    if [ ! -s $input/$snv_file ]
    then
        $script_path/errorlog.sh $input/$snv_file polyphen.sh ERROR "not found"
		exit 1;
    fi
    
    cat $input/$snv_file | awk '$0 !~ /^#/' | awk '{print $1":"$2"\t"$4"/"$5}' > $polyphen/$snv_file.poly
    num=`cat $polyphen/$snv_file.poly |wc -l `
    ### get the uniport id
    if [ $num -gt 0 ]
    then
        $script_path/parallel.poly.pl $threads $polyphen/$snv_file.poly $genome_build $pph $polyphen/$snv_file.poly.uniprot
    else
        touch $polyphen/$snv_file.poly.uniprot
    fi
    rm $polyphen/$snv_file.poly
    
    num=`cat $polyphen/$snv_file.poly.uniprot | wc -l`
    if [ $num -gt 0 ]
    then
        ### get the prediction
        uniprot=`awk '{ for(i=1;i<=NF;i++){ if ($i == "spacc") {print i} } }' $polyphen/$snv_file.poly.uniprot`
        aa1=`awk '{ for(i=1;i<=NF;i++){ if ($i == "aa1") {print i} } }' $polyphen/$snv_file.poly.uniprot`
        aa2=`awk '{ for(i=1;i<=NF;i++){ if ($i == "aa2") {print i} } }' $polyphen/$snv_file.poly.uniprot`
        aapos=`awk '{ for(i=1;i<=NF;i++){ if ($i == "cdnpos") {print i} } }' $polyphen/$snv_file.poly.uniprot`
        cat $polyphen/$snv_file.poly.uniprot | awk 'NR>1' | awk -v a1=$aa1 -v a2=$aa2  -F '\t' '$a1 !~ $a2' |  awk -v uni=$uniprot 'length($uni)>1' | awk -v chr=$snp_pos -v uni=$uniprot -v a1=$aa1 -v a2=$aa2 -v pos=$aapos -F '\t' '{if ($a1 !~ "*" && $a2 !~ "*") print $uni"\t"$pos"\t"$a1"\t"$a2};' | awk -F '\t' 'NF==4' > $polyphen/$snv_file.poly.uniprot.in
        $script_path/split.pl $threads $polyphen/$snv_file.poly.uniprot.in $genome_build $pph $polyphen/$snv_file.poly.uniprot.in.predict
        rm $polyphen/$snv_file.poly.uniprot.in
    else
        touch $polyphen/$snv_file.poly.uniprot.in.predict
    fi    
    $script_path/map.polyphen.pl $polyphen/$snv_file.poly.uniprot $polyphen/$snv_file.poly.uniprot.in.predict > $polyphen/$snv_file.poly
    rm $polyphen/$snv_file.poly.uniprot
    rm $polyphen/$snv_file.poly.uniprot.in.predict
    if [ ! -f $polyphen/$snv_file.poly ]
	then
		$script_path/errorlog.sh $polyphen/$snv_file.poly polyphen.sh ERROR "not created"
		exit 1;
	fi	
	echo `date`
fi    
    
    
    
