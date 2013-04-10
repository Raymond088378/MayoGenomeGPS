#!/bin/bash

if [ $# != 1 ]
then
	echo -e "script to check the config files\nUsage:  <full path to run info file>"
	exit 1;
fi	
	
run_info=$1
dir=`dirname $run_info`
touch $dir/config.log 
if [ ! -s $run_info ]
then
	echo -e "ERROR : run_info=$run_info does not exist.\n" >> $dir/config.log
	exit 1;
else
	dos2unix $run_info
	## removing trailing and leading spaces from run info file
	cat $run_info | sed 's/^[ \t]*//;s/[ \t]*$//' > $run_info.tmp
	mv $run_info.tmp $run_info
fi


### check for tool info file
tool_info=$( cat $run_info | grep -w '^TOOL_INFO' | cut -d '=' -f2)
if [ ! -s $tool_info ]
then
	echo -e "ERROR : tool_info=$tool_info does not exist \n" >> $dir/config.log
	exit 1;
else
	dos2unix $tool_info
	cat $tool_info | sed 's/^[ \t]*//;s/[ \t]*$//' > $tool_info.tmp
	mv $tool_info.tmp $tool_info
fi

### check for sample info file
sample_info=$( cat $run_info | grep -w '^SAMPLE_INFO' | cut -d '=' -f2)
if [ ! -s $sample_info ]
then
	echo -e "ERROR : sample_info=$sample_info does not exist \n" >> $dir/config.log
	exit 1;
else
	dos2unix $sample_info
	cat $sample_info | sed 's/^[ \t]*//;s/[ \t]*$//' > $sample_info.tmp
	mv $sample_info.tmp $sample_info
fi

### check for memory info file
memory_info=$( cat $run_info | grep -w '^MEMORY_INFO' | cut -d '=' -f2)
if [ ! -s $memory_info ]
then
	echo -e "ERROR : memory_info=$memory_info does not exist \n" >> $dir/config.log
	exit 1;
else
	dos2unix $memory_info
	cat $memory_info | sed 's/^[ \t]*//;s/[ \t]*$//' > $memory_info.tmp
	mv $memory_info.tmp $memory_info
fi
