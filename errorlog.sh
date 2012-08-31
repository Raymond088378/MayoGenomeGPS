#!/bin/bash

if [ $# != 4 ]
then
	echo -e "srcipt to report error in recommnended format \n <file name><script_name><ERROR/WARNING> <free tesxt>"
else
	file=$1
	script=$2
	type=`echo $3 | tr "[a-z]" "[A-Z]"`
	what=$4   ## empty or missing or free text
	echo "$type : $file file is $what in script $script [`date`]"
fi
