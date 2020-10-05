#!/bin/bash

scriptname=`basename $0`
tmpfile=$(mktemp ./${basename}.$$)

cat .github/workflows/python-package.yml | grep black  | grep check | sed -e "s/^\s+//" -e "s/\"/'/g"  >  $tmpfile
cat .github/workflows/python-package.yml | grep flake8 | grep count | sed -e "s/^\s+//" -e "s/10/200/" >> $tmpfile

echo "Run these commands:"
cat $tmpfile

while IFS= read command
do
  echo "Command 1: #${command}#"
  command=`echo $command | sed -e "s/\s+//"`
  echo "Command 2: #${command}#"
  $command
done < "$tmpfile"

rm -f "$tmpfile"

exit 0
