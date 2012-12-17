#!/bin/bash

if [ $# != 1 ]
then
	echo -e "script to wait for the fix \nusage: ./wait.sh </path/to/filename>"
else
	echo `date`
	file=$1
	while [ -f $file ]
	do
		sleep 5m
		echo "waiting for the $file to be fixed"
	done
	echo `date`
fi			