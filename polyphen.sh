#!/bin/sh

if [ $# != 4 ]
then
    echo "Usage: wrapper to run polyphen \n <polphen dir> <input directory><sample><run info>"
else
    set -x
    echo `date`
    polyphen=$1
    input=$2
    sample=$3
    run_info=$4
    
    chr=$(cat $run_info | grep -w '^CHRINDEX' | cut -d '=' -f2 | tr ":" "\n" | head -n $SGE_TASK_ID | tail -n 1)
    tool_info=$(cat $run_info | grep -w '^TOOL_INFO' |  cut -d '=' -f2)
    pph=$(cat $tool_info | grep -w '^POLYPHEN' |  cut -d '=' -f2)
    perl_lib=$(cat $tool_info | grep -w '^PERL_POLYPHEN_LIB' |  cut -d '=' -f2)
    genome_build=$(cat $run_info | grep -w '^GENOME_BUILD' |  cut -d '=' -f2)
    perl=$(cat $tool_info | grep -w '^PERL_CIRCOS' |  cut -d '=' -f2)
    
    export PERL5LIB=${PERL5LIB}:$perl_lib
    
    snv_file=$sample.chr$chr.SNV.filter.i.c.vcf
    
    cat $input/$snv_file | awk '$0 !~ /^#/' | awk '{print $1":"$2"\t"$4"/"$5}' > $polyphen/$snv_file.poly
    
    ### get the uniport id
    $perl $pph/bin/mapsnps.pl -g $genome_build -v 0 $polyphen/$snv_file.poly > $polyphen/$snv_file.poly.uniprot
    rm $polyphen/$snv_file.poly
    ### get the prediction
    uniprot=`awk '{ for(i=1;i<=NF;i++){ if ($i == "spacc") {print i} } }' $polyphen/$snv_file.poly.uniprot`
    aa1=`awk '{ for(i=1;i<=NF;i++){ if ($i == "aa1") {print i} } }' $polyphen/$snv_file.poly.uniprot`
    aa2=`awk '{ for(i=1;i<=NF;i++){ if ($i == "aa2") {print i} } }' $polyphen/$snv_file.poly.uniprot`
    aapos=`awk '{ for(i=1;i<=NF;i++){ if ($i == "cdnpos") {print i} } }' $polyphen/$snv_file.poly.uniprot`
    
    cat $polyphen/$snv_file.poly.uniprot | awk 'NR>1' | awk -v a1=$aa1 -v a2=$aa2  -F '\t' '$a1 !~ $a2' |  awk -v uni=$uniprot 'length($uni)>1' | awk -v chr=$snp_pos -v uni=$uniprot -v a1=$aa1 -v a2=$aa2  \
        -v pos=$aapos -F '\t' '{print $uni"\t"$pos"\t"$a1"\t"$a2}' > $polyphen/$snv_file.poly.uniprot.in
    
    $perl $pph/bin/run_pph.pl -g $genome_build -d $polyphen -v 0 $polyphen/$snv_file.poly.uniprot.in > $polyphen/$snv_file.poly.uniprot.in.predict
    rm $polyphen/$snv_file.poly.uniprot.in
    
    ### map chr pos with the prediction
    
    perl $script_path/map.polyphen.pl $polyphen/$snv_file.poly.uniprot $polyphen/$snv_file.poly.uniprot.in.predict > $polyphen/$snv_file.poly
    rm $polyphen/$snv_file.poly.uniprot
    rm $polyphen/$snv_file.poly.uniprot.in.predict
    echo `date`
fi    
    
    
    