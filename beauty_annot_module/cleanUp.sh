#! /bin/bash
## Remove directory only if it begins with "tmp" and is followed immediately by a digit.

DIR2REMOVE=$(basename $1)
echo "attempt to remove: $DIR2REMOVE"
if [[ $DIR2REMOVE =~ ^tmp[0-9] ]] ; then
    echo "OK"
	rm -R $1
else
    echo "not OK"
fi