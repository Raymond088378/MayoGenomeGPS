#!/bin/sh

if [ $# != 4 ]
then
	echo "<file name><script_name>"
else	
	file=$1
	script=$2
	type=`$3 | tr "[a-z]" "[A-Z]"`
	what=$4   ## empty or missing
	echo " $type : $file is $what in script $script [`date`]"
fi	

