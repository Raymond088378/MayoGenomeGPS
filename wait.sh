#!/bin/bash

if [ $# != 1 ]
then
	echo -e "script to wait for the fix \nusage: <filename>"
else
	file=$1
	while [ -f $file ]
	do
		sleep 2m
		echo "waiting for the $file to be fixed"
	done
fi			