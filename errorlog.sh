#!/bin/bash

if [ $# != 4 ]
then
	echo -e "script to report error in recommnended format\nUsage: errorlog.sh <file name><script_name><ERROR/WARNING> <free tesxt>"
else
	file=$1
	script=$2
	type=`echo $3 | tr "[a-z]" "[A-Z]"`
	what=$4   ## empty or missing or free text
	echo "$type : $file file is $what in script $script [`date`]"
fi
