#!/bin/bash

if [ $# != 4 ];
then
	echo -e "script to add allele frequency for hapmap and 1kgenome populations\nUsage ./add_3population.sh [input file with chr as 1st column and position as 2nd column, 1-based] [outout name, e.g., run70 ]"
else
	set -x
	echo `date`
	run_info=$3
	chr=$4
	tool_info=$( cat $run_info | grep -w '^TOOL_INFO' | cut -d '=' -f2)
	script_path=$( cat $tool_info | grep -w '^WORKFLOW_PATH' | cut -d '=' -f2)
	hapmap=$( cat $tool_info | grep -w '^HAPMAP' | cut -d '=' -f2)
	kgenome=$( cat $tool_info | grep -w '^KGENOME' | cut -d '=' -f2)
	GenomeBuild=$( cat $run_info | grep -w '^GENOMEBUILD' | cut -d '=' -f2)
	cat $hapmap/all_allele_freqs_CEU.txt  | grep -w chr$chr > $1.all_allele_freqs_CEU.txt
	cat $hapmap/all_allele_freqs_YRI.txt  | grep -w chr$chr > $1.all_allele_freqs_YRI.txt
	cat $hapmap/all_allele_freqs_JPT+CHB.txt  | grep -w chr$chr > $1.all_allele_freqs_JPT+CHB.txt
	cat $kgenome/CEU.$GenomeBuild | grep -w chr$chr > $1.CEU.$GenomeBuild
	cat $kgenome/YRI.$GenomeBuild | grep -w chr$chr > $1.YRI.$GenomeBuild
	cat $kgenome/JPT+CHB.$GenomeBuild | grep -w chr$chr > $1.JPT+CHB.$GenomeBuild
	
	perl $script_path/add_hapmap_1kgenome_allele_freq.pl -i $1 -c 1 -p 2 -b 1 -r $chr -e CEU -s $1.all_allele_freqs_CEU.txt -g $1.CEU.$GenomeBuild -o $2.CEU&&$perl $script_path/add_hapmap_1kgenome_allele_freq.pl -i $2.CEU -c 1 -r $chr -p 2 -b 1 -e YRI -s $1.all_allele_freqs_YRI.txt -g $1.YRI.$GenomeBuild -o $2.CEU.YRI&&perl $script_path/add_hapmap_1kgenome_allele_freq.pl -i $2.CEU.YRI -c 1 -p 2 -r $chr -b 1 -e JPT+CHB -s $1.all_allele_freqs_JPT+CHB.txt -g $1.JPT+CHB.$GenomeBuild -o $2.CEU.YRI.CHBJPT.txt
	rm $1.all_allele_freqs_CEU.txt $1.all_allele_freqs_YRI.txt $1.all_allele_freqs_JPT+CHB.txt $1.CEU.$GenomeBuild $1.YRI.$GenomeBuild $1.JPT+CHB.$GenomeBuild
	echo `date`
fi